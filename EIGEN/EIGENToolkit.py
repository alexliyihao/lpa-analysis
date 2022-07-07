import pandas as pd
import textwrap
import os


class EIGENToolkit():
    """A class generate EIGENSTRAT format input for EIGENSOFT

    This class intake DataFrame from allele encoding table and provides the following APIs:

    Initialization:
        etk = EIGENToolkit(chromosome = 6,
                           genetic_position = 0,
                           physical_position_base = 161033785,
                           verbose = 1)
        chromosome: int, the number of chromosome used for encoding
        genetic_position, physical_position_base: check EIGENSOFT's documentation
        verbose: 1 or 0, if 1, when new file is generated a notice will be printed
    One-line Output:
        etk.to_eigen(snp, ethnicity_from_pca, gender, output_path, filename)
            the following three pandas instance will be inner-join-aligned:
                snp: pd.DataFrame, the snp integer-encoding table,
                     mutations as columns and individual sample as rows
                ethnicity_from_pca: pd.DataFrame, the ethnicity determined by eigenstrat
                gender: pd.Series, the gender of each individuals
            output_path: str, the output path, if not exist, it will create one
            filename: str, the file name for all the eigenstrat input without extension name
    Reading helper:
        etk.read_PC_weight(path): read a .snpweight output file
        etk.read_PCA_result(path): read a .pca.evec output file
    """

    def __init__(self,
                 chromosome = 6,
                 genetic_position = 0,
                 physical_position_base = 161033785,
                 verbose = 1):
        self._chromosome = chromosome
        self._genetic_position = genetic_position
        self._physical_position_base = physical_position_base
        self._verbose = verbose

    def _snp_pos_to_eigenstrat(self, snp_pos):
        """
        convert a snp_pos format e.g."11-G/A" to eigenstrat format .snp row
        """
        # split position and ref/alt
        pos, snp = snp_pos.split("-")
        # split ref and alt
        ref, alt = snp.split("/")
        # the actual pos start at 161033785 as 1, so
        actual_pos = self._physical_position_base + int(pos) - 1
        output = [snp_pos, self._chromosome, self._genetic_position, actual_pos, ref, alt]
        return "\t".join([str(i) for i in output])

    def _generate_eigenstrat_snp(self, snp, output):
        """
        generate eigenstrat .snp file from a snp pandas.DataFrame
        """
        snp_file = [snp_pos_to_eigenstrat(snp_pos = snp_pos) for snp_pos in snp.columns]
        with open(output, "w") as file:
            file.write("\n".join(snp_file))

    def _generate_eigenstrat_ind(self, complete_table, output):
        """
        generate a ind file for individual information
        """
        complete_table = complete_table[["GENDER", "ethnicity"]].copy()
        complete_table.loc[:,"GENDER"] = complete_table["GENDER"].apply(lambda x: "M" if x else "F")
        complete_table.to_csv(output, sep='\t', header=False)

    def _generate_eigenstrat_geno(self, snp, output):
        """
        generate a .eigenstratgeno file for eigenstrat
        """
        snp = 2-snp.T
        with open(output, "w") as file:
            file.writelines(
                "".join(str(i) for i in snp.iloc[i,:].values)+"\n"
                for i in range(snp.shape[0])
            )

    def _generate_eigenstrat_param(self, filename):
        """
        generate a string for .pca.par file for eigenstrat
        """
        return textwrap.dedent(f"""\
        genotypename: {filename}.eigenstratgeno
        snpname: {filename}.snp
        indivname: {filename}.ind
        evecoutname: {filename}.pca.evec
        evaloutname: {filename}.eval
        altnormstyle: NO
        numoutevec: 3
        numoutlieriter: 5
        numoutlierevec: 10
        outliersigmathresh: 6
        qtmode: 0
        snpweightoutname: {filename}.snpweight
        """)

    def to_eigen(snp, ethnicity_from_pca, gender, output_path, filename):
        """
        One-line wrapper all the previous code, the output should be a path rather than a file
        Args:
            snp: pd.DataFrame, the snp table
            ethnicity_from_pca: pd.DataFrame, the ethnicity determined by eigenstrat
            gender: pd.Series, the gender of each individuals
            output_path: str, the output path, if not exist, will create one
            filename: str, the file name for all the eigenstrat input without extension name
        Saves:
            {filename}.snp,
            {filename}.ind,
            {filename}.eigenstratgeno,
            {filename}.pca.par
        """
        # create the output path
        os.makedirs(output_path, exist_ok = True)

        snp,ethnicity_from_pca = snp.align(ethnicity_from_pca, axis=0, join = "inner")
        snp,gender = snp.align(gender, axis=0, join = "inner")
        complete_table = pd.concat([snp, ethnicity_from_pca, gender], axis = 1)
        print(f"complete dataset: sample size {snp.shape[0]}, snp pool size {snp.shape[1]}")
        self._generate_eigenstrat_snp(
            snp = snp,
            output = os.path.join(output_path, f"{filename}.snp")
        )
        if self._verbose == 1:
            print(f"{filename}.snp saved")

        self._generate_eigenstrat_ind(
            complete_table = complete_table,
            output = os.path.join(output_path, f"{filename}.ind")
        )
        if self._verbose == 1:
            print(f"{filename}.ind saved")

        self._generate_eigenstrat_geno(
            snp = snp,
            output = os.path.join(output_path, f"{filename}.eigenstratgeno")
        )
        if self._verbose == 1:
            print(f"{filename}.eigenstratgeno saved")

        param_string = self._generate_eigenstrat_param(filename)
        with open(os.path.join(output_path, f"{filename}.pca.par"), "w") as file:
            file.write(param_string)
        if self._verbose == 1:
            print(f"{filename}.pca.par saved")

    def read_PC_weight(self, path):
        """read a .snpweight output file and tidy the format"""
        eigenstrat_result = pd.read_csv(full_path,
                                        sep = "\s+",
                                        header = None)
        eigenstrat_result = eigenstrat_result.drop(columns = [1,2])
        eigenstrat_result = eigenstrat_result.rename(columns = {0: "ID",
                                                                3: "PC1_weight",
                                                                4: "PC2_weight",
                                                                5: "PC3_weight"})
        return eigenstrat_result

    def read_PCA_result(self, path):
        """read a .pca.evec output file and tidy the format"""
        eigenstrat_result = pd.read_csv(full_path,
                                 skiprows = 1,
                                 sep = "\s+",
                                 header = None)
        eigenstrat_result = eigenstrat_result.rename(columns = {0: "ID",
                                                                1: "PC1",
                                                                2: "PC2",
                                                                3: "PC3",
                                                                4: "ethnicity"})
        eigenstrat_result["ethnicity"] = eigenstrat_result[["ethnicity"]].replace([1,2,3], ["EU", "AF", "HISP"])
        return eigenstrat_result
