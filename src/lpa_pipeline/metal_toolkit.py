"""The toolkit running METAL-related pipeline internally

Example:

    The class should be initiated as follows::

        mtk = metal_toolkit.METALToolkit(
            ethnicity = ["EU", "AF", "HISP"],
            engine_list = ["Logit", "OLS"],
            verbose = 1,
            metal_path = "/mnt/mfs/cluster/bin/METAL/metal",
            snp_alias = "variant"
            )

    One-line Pipeline::

        mtk.run_metal(
            path_association = "/path/of/association/module/output",
            path_meta = "/path/of/meta_analysis/output")

    It internally does everything below:

        mtk.copy_tree(path_association = path_association, path_meta = path_meta)
            move all the file in path_association to path_meta,
            and properly arrange them by traits

        mtk.csv_to_metal(path)
            Prepare a METAL format CSV file from existing CSV file at <path>

        mtk.generate_metal_script(path)
            Prepare  metal.sh and metal_script.txt with specified ethnicity.

        subprocess("chmod +x metal.sh" & "metal source metal_script.txt")
            run METAL from the script generated

        mtk.aggregate_result(path):
            After running the METAL analysis, this API will generate an
            aggregate table for this folder for next step analysis

        pd.concat([path/aggregate.csv for path in sub_folders], axis = 0)
            concatenate all the aggregate.csv

        mtk.formatting_final_result(df, multi_line_header = False)
            format the aggregate.csv
            multi_line_header is a historical problem, use False is fine.
            Set to True will generate another header row.

        mtk.save_final_output(df,path = "path", multi_line_header = False)
            deals with a bug in pandas when using multi-line headers

The ``multi_line_header`` only makes difference at the output header,
it will create another header row, just for the visualization pipeline
"""

import pandas as pd
import os
import shutil
import glob
import subprocess
import shlex
import re
from typing import List


class METALToolkit:
    """Pipeline for python preparing/processing METAL's meta analysis"""

    def __init__(self,
                 ethnicity: List[str] = ["EU", "AF", "HISP"],
                 engine_list: List[str] = ["OLS", "Logit"],
                 verbose: int = 1,
                 metal_path: str = "/mnt/mfs/cluster/bin/METAL/metal",
                 snp_alias: str = "variant") -> None:
        """initializer save the ethnicity setting,

        Args:
            ethnicity: list of str, the list of ethnicities in string format,
                        if not provided, using ["eu", "af", "hisp"]
            engine_list: list of str, the class name of all the engine included, used for filename processing
            verbose: int, the verbosity running the pipeline.
            metal_path: str, the path to METAL binary, the development environment cannot use
            "module load METAL" via Jupyterlab hold in SGE, you can either
            1. "module load METAL && echo $PATH" find the path
            2. override generate_metal_running_script() method when your environment allows.
            snp_alias: str, the setting from association module
        """
        self._ethnicity = [i.replace(" ", "") for i in ethnicity]
        self._verbose = verbose
        self._metal_app_path = metal_path
        self._snp_alias = snp_alias
        self._engine_list = engine_list

    def run_metal(
            self,
            path_association: str,
            path_meta: str,
            rename: bool = False,
            multi_line_header: bool = False
    ) -> pd.DataFrame:
        """API processing all the tables at once in a folder

        It will copy the path from <path_association> (keep an unaffected
        version) to <path_meta> (must not already exist!) and process
        everything in new path to

        Args:
            path_association: str, the path of association outputs
            path_meta: str, the path to meta-analysis, it must not already exist.
            rename: boolean, using False is fine, just for formatting
            multi_line_header: boolean, using False is fine.

        Returns:
            pd.DataFrame: an aggregate result for meta-analysis

        Saves:
            in path_meta:
            - csv files with an additional "_for_metal" at the end of file name
            - metal.sh and metal_script.txt for metal running
            - aggregate.csv in each trait folder
            - aggregate_all_traits.xlsx in path_meta for better exploration.
        """
        self._logging("copying the association file...")
        self.copy_tree(path_association=path_association,
                       path_meta=path_meta)

        self._logging("Generating the association file for METAL...")
        file_list = glob.glob(os.path.join(
            path_meta, "**", "*.csv"), recursive=True)
        for tables in tqdm(file_list):
            self.csv_to_metal(tables)

        self._logging("Generating the METAL scripts...")
        self.generate_metal_script(path=path_meta)

        self._logging("Running METAL...")
        subprocess.check_call(shlex.split(
            f"chmod +x {os.path.join(path_meta, 'metal.sh')}"))
        subprocess.check_output(shlex.split(
            f"sh {os.path.join(path_meta, 'metal.sh')}"))

        self._logging("Aggregating the output in each ethnicity...")
        sub_folders = glob.glob(os.path.join(path_meta, "**", ""))
        for trait_folder in tqdm(sub_folders):
            aggregate_df = self.aggregate_result(path=trait_folder)
            aggregate_df.to_csv(os.path.join(trait_folder, "aggregate.csv"))

        self._logging("Aggregating the output from each traits...")
        all_result = pd.concat((pd.read_csv(os.path.join(i, "aggregate.csv"))
                                for i in tqdm(sub_folders)),
                               axis=0)
        all_result = self.formatting_final_result(
            all_result,
            multi_line_header=multi_line_header,
            rename=rename)
        self.save_final_output(
            df=all_result,
            path=os.path.join(path_meta, "aggregate_all_traits.xlsx"),
            multi_line_header=multi_line_header)
        return all_result

    # --------------------------prepare the metal table-----------------------------

    def dataframe_to_metal(self, df: pd.DataFrame) -> pd.DataFrame:
        """for a given association DataFrame, prepare a METAL accepted version

        Args:
            df: pd.DataFrame, an association table generated by
            association.SNPAssociation instance

        Returns:
            pd.DataFrame: a tidied version of association table
            ready for METAL's meta analysis
        """
        df_new = pd.DataFrame()
        df_new["SNP"] = df[self._snp_alias].apply(self._get_pos)
        df_new["REF_ALLELE"] = df[self._snp_alias].apply(
            lambda i: i.split("-")[-1]).apply(lambda i: i.split("/")[0])
        df_new["ALT_ALLELE"] = df[self._snp_alias].apply(
            lambda i: i.split("-")[-1]).apply(lambda i: i.split("/")[-1])
        df_new["FREQ"] = df["rel_freqs"]
        df_new["BETA"] = df[f"{self._snp_alias} Beta"]
        df_new["SE"] = df[f"{self._snp_alias} SE"]
        df_new["PVAL"] = df[f"{self._snp_alias} P-value"]
        df_new["TOTAL"] = df["n_sample"]
        df_new["COUNT"] = df["abs_freqs"]
        df_new["POS-REF/ALT"] = df[self._snp_alias]
        return df_new

    def csv_to_metal(self, path: str):
        """API including the csv file reading and writing

        Args:
            path: string, the path of an association table generated by
            ...pipeline.association.SNPAssociation instance
        Saves:
            csv file with an additional "_for_metal" at the end of the file name
        """
        df = pd.read_csv(path).rename(columns={"Unnamed: 0": self._snp_alias})
        df_new = self.dataframe_to_metal(df)
        new_path = path.replace(".csv", "_for_metal.csv")
        df_new.to_csv(new_path, index=False)

    # --------------------------prepare the metal script----------------------------

    def generate_metal_script(self, path=None) -> None:
        """API generate the metal script for the meta analysis,

        Args:
            path: Optional[string], if given, it will save the script generated
            to the path given as well
        Returns:
            str: the script for running the meta analysis
        """
        metal_script_text = "\n".join([self._metal_script_header(),
                                       self._metal_script_description(),
                                       self._metal_script_analyze()])
        running_script_text = self.generate_metal_running_script(
            metal_path=path)
        if path is not None:
            with open(os.path.join(path,
                                   "metal_script.txt"), "w") as metal_script:
                metal_script.write(metal_script_text)
            with open(os.path.join(path, "metal.sh"), "w") as running_script:
                running_script.write(running_script_text)
            if self._verbose == 1:
                print(f"METAL scripts are saved to {path}")

    def _metal_script_header(self) -> str:
        """prepare a header for the metal script"""
        return textwrap.dedent("""\
        SEPARATOR  COMMA
        # Meta-analysis weighted by standard error does not work well
        # when different studies used very different transformations.
        # In this case, some attempt was made to use similar trait
        # transformation and you can request a standard error based
        # analysis by uncommenting the following line:
        SCHEME   STDERR

        # Usually, it is a good to apply genomic control to each
        # input file. However, in this example, all the markers being
        # examined reside in strongly associated loci and it probably
        # is not a good idea. To find out what happens if you turn
        # on genomic control, uncomment the following line.
        # GENOMICCONTROL ON

        # To help identify allele flips, it can be useful to track
        # allele frequencies in the meta-analysis. To enable this
        # capability, uncomment the following two lines.
        # AVERAGEFREQ ON
        # MINMAXFREQ ON

        # To restrict meta-analysis to two previously reported SNPs
        # and summarize study specific results, uncomment the two
        # lines that follow.
        # ADDFILTER SNP IN (rs10830963,rs563694)
        # VERBOSE ON

        """)

    def _metal_script_description(self) -> str:
        """prepare the description part of METAL script"""
        return "\n".join([textwrap.dedent(f"""\
            # Describe and process the DGI input files

            MARKER   POS-REF/ALT
            ALLELE   ALT_ALLELE REF_ALLELE
            FREQ     FREQ
            EFFECT   BETA
            STDERR   SE
            PVAL     PVAL

            PROCESS ethnicity={ethnicity}_for_metal.csv
            """) for ethnicity in self._ethnicity])

    def _metal_script_analyze(self):
        """prepare the analysis part of METAL script"""
        return textwrap.dedent(f"""\
            # Execute meta-analysis
            ANALYZE
            """)

    def generate_metal_running_script(self, metal_path) -> str:
        """prepare the running script of METAL"""
        return textwrap.dedent(
            f"""\
            for trait in {os.path.join(metal_path, "*", "")} ; do
                cd "$trait";
                {self._metal_app_path} source {os.path.join(metal_path, "metal_script.txt")} > metal_$(basename "$trait").log;
                cd {metal_path};
            done
            """)

    # ---------------------aggregate the metal result-------------------------

    def _get_snp_data(self, path: str, key: str) -> pd.DataFrame:
        """for a folder, and a specific key, get the requested info, add suffix"""
        df = pd.read_csv(os.path.join(
            path, f"ethnicity={key}.csv"), index_col=0)
        df = df[[f"{self._snp_alias} Beta", f"{self._snp_alias} SE",
                 f"{self._snp_alias} P-value",
                 "rel_freqs", "abs_freqs", "n_sample"]]
        df.columns = df.columns + f"_{key}"
        return df

    def _get_pos(self, variant: str) -> int:
        """extract pos in integer from a string format pos-ref/alt"""
        return int(variant.split("-")[0])

    def aggregate_result(self, path: str) -> pd.DataFrame:
        """give a folder with METAL result, combine all the results"""
        meta = pd.read_csv(os.path.join(
            path, "METAANALYSIS1.TBL"), sep="\t", index_col=0)
        test_list = [meta] + [self._get_snp_data(path, key) for key in self._ethnicity]
        final = pd.concat(test_list, axis=1, join="outer")
        final = final.reset_index()
        final["pos"] = final["index"].apply(self._get_pos)
        final = final.sort_values("pos").set_index("index").drop(columns="pos")
        final["trait"] = path.split("/")[-2]
        return final

    def formatting_final_result(
            self,
            df: pd.DataFrame,
            multi_line_header: bool = False,
            rename: bool = False
    ) -> pd.DataFrame:
        df["total_count"] = sum((df[f"abs_freqs_{i}"].fillna(0) for i in self._ethnicity))
        df["total_population"] = sum((df[f"n_sample_{i}"].fillna(0) for i in self._ethnicity))
        colnames = [[f"{self._snp_alias} Beta_{i}", f"{self._snp_alias} SE_{i}",
                     f"{self._snp_alias} P-value_{i}",
                     f"rel_freqs_{i}", f"abs_freqs_{i}", f"n_sample_{i}"] for i in self._ethnicity]
        df = df[["index"] +
                [item for sublist in colnames for item in sublist] +
                ["Allele1", "Allele2", "Effect", "StdErr", "P-value",
                 "Direction", "total_count", "total_population", "trait"]
                ]
        if rename:
            df.rename(columns={
                "index": ("SNPs", "MarkerName"),
                f"{self._snp_alias} Beta_EU": ("EUR", f"{self._snp_alias}.Beta.EUR"),
                f"{self._snp_alias} SE_EU": ("EUR", f"{self._snp_alias}.SE.EUR"),
                f"{self._snp_alias} P-value_EU": ("EUR", f"{self._snp_alias}.P.value.EUR"),
                "rel_freqs_EU": ("EUR", "rel_freqs.EUR"),
                "abs_freqs_EU": ("EUR", "abs_freqs.EUR"),
                "n_sample_EU": ("EUR", "n_sample.EUR"),
                f"{self._snp_alias} Beta_AF": ("AFR", f"{self._snp_alias}.Beta.AA"),
                f"{self._snp_alias} SE_AF": ("AFR", f"{self._snp_alias}.SE.AA"),
                f"{self._snp_alias} P-value_AF": ("AFR", f"{self._snp_alias}.P.value.AA"),
                "rel_freqs_AF": ("AFR", "rel_freqs.AA"),
                "abs_freqs_AF": ("AFR", "abs_freqs.AA"),
                "n_sample_AF": ("AFR", "n_sample.AA"),
                f"{self._snp_alias} Beta_HISP": ("HISP", f"{self._snp_alias}.Beta.HISP"),
                f"{self._snp_alias} SE_HISP": ("HISP", f"{self._snp_alias}.SE.HISP"),
                f"{self._snp_alias} P-value_HISP": ("HISP", f"{self._snp_alias}.P.value.HISP"),
                "rel_freqs_HISP": ("HISP", "rel_freqs.HISP"),
                "abs_freqs_HISP": ("HISP", "abs_freqs.HISP"),
                "n_sample_HISP": ("HISP", "n_sample.HISP"),
                "Allele1": ("META", "Allele1"),
                "Allele2": ("META", "Allele2"),
                "Effect": ("META", "Effect"),
                "StdErr": ("META", "StdErr"),
                "P-value": ("META", "P.value"),
                "Direction": ("META", "Direction"),
                "total_count": ("META", "META.ALLELE.COUNT"),
                "total_population": ("META", "META.ALLELE.N"),
                "trait": ("META", "trait")
            },
                inplace=True)
            df.columns = pd.MultiIndex.from_tuples(df.columns)
            if multi_line_header is False:
                df.columns = df.columns.droplevel(0)
        df.reset_index(drop=True, inplace=True)
        return df

    # --------------------------File Management---------------------------------
    def copy_tree(self, path_association: str, path_meta: str) -> None:
        """move file from path_association to path_meta, arrange by traits"""
        # copy the path from <path_association> to <path_meta>
        shutil.rmtree(path_meta, ignore_errors=True)
        shutil.copytree(src=path_association, dst=path_meta)
        for filename in tqdm(os.listdir(path_meta)):
            info = filename.split("_")
            if info[2] == "OLS":
                trait = f"{info[0]}_{info[1]}"
                ethnicity = info[3]
            else:
                trait = info[0]
                ethnicity = info[2]
            for engine_name in self._engine_list:
                if engine_name in filename:
                    trait = re.search(f"^.*?(?=_{engine_name}_)", filename).group(0)
                    ethnicity = re.search(f"(?<=_{engine_name}_)(.*)(?=_N_snp)", filename).group(0)
            os.makedirs(os.path.join(path_meta, trait), exist_ok=True)
            if os.path.splitext(filename)[1] == ".csv":
                shutil.move(os.path.join(path_meta, filename),
                            os.path.join(path_meta, trait, f"{ethnicity}.csv").replace(" ", ""))

    def save_final_output(
            self,
            df: pd.DataFrame,
            path: str,
            multi_line_header: bool
    ) -> None:
        """
        wrapper saving the output file, which is for solving a bug saving
        pd.DataFrame with multiIndex to excel,
        see https://stackoverflow.com/a/71305025.
        """
        if multi_line_header is True:
            with pd.ExcelWriter(path, engine="xlsxwriter") as xl_writer:
                self._save_double_column_df(df=df, xl_writer=xl_writer)
        else:
            df.to_excel(path)

    def _logging(self, log: str) -> None:
        """just for code cleaning..."""
        if self._verbose is True:
            print(log)

    def _save_double_column_df(
            self,
            df: pd.DataFrame,
            xl_writer: pd.ExcelWriter,
            start_row: int = 0,
            **kwargs) -> None:
        """Function saving two-level column DataFrame to xlwriter
            credit to https://stackoverflow.com/a/71305025

        Args:
             df: pd.DataFrame, the table to save
             xl_writer: pd.ExcelWriter, book for saving
             start_row: int, row from which dataframe will begin
             **kwargs: Any, arguments of `to_excel` function of DataFrame`
        """
        # inputs:

        df.drop(df.index).to_excel(xl_writer, startrow=start_row, **kwargs)
        df.to_excel(xl_writer, startrow=start_row + 1, header=False, **kwargs)


if __name__ == "__main__":
    import argparse
    import textwrap
    parser = argparse.ArgumentParser(prog="metal_toolkit.py",
                                     description=textwrap.dedent("""
        This script copies the result of Association.py to a new folder
          - The purpose for this copying is to keep a untouched version of data.
          - Run METAL internally, all the scripts are generated by the code.
          - Aggregate all the results
          - The output is a .xlsx file rather than a csv file for exploratory purpose
        """), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-I", "--input_path", type=str, required=True,
                        help="path to the folder storing association pipeline output")
    parser.add_argument("-O", "--output_path", type=str, required=True,
                        help="path saving the output and intermediate results, this path must not already exist")
    parser.add_argument("-L", "--ethnicity_list", required=True, nargs="+",
                        help="the ethnicity mentioned, case sensitive, used as '-L EU AF HISP'")
    parser.add_argument("-P", "--metal_path", type=str, required=True,
                        help="""the path to METAL binary, you can either
                1. "module load METAL && echo $PATH" find the path
                2. modify generate_metal_running_script() method.""")
    parser.add_argument("-V", "--verbosity", type=int, choices=[0, 1],
                        required=False, default=1,
                        help="verbosity 0 or 1, if set to 1, will print some logs, default 1")
    parser.add_argument("-M", "--multi_line_header", required=False,
                        action="store_true",
                        help="If mentioned, the output will have another header row, just for formatting")
    Args = parser.parse_args()
    mtk = METALToolkit(ethnicity=Args.ethnicity,
                       verbose=Args.verbosity,
                       metal_path=Args.metal_path)
    mtk.run_metal(path_association=Args.input_path,
                  path_meta=Args.output_path,
                  multi_line_header=Args.multi_line_header)
