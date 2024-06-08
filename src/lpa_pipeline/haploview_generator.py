"""A class generating output for Haploview4.2 linkage format

Common usage:

Given two pandas.DataFrame

* ``encoded_result`` as the output of ``snps_filter`` module
* ``locus_table`` as the output of ``locus_collector`` module
* ``ethnicity_info`` as a pd.DataFrame with a column ``ethnicity`` including the ethnicity info for each individual

    Initialize::

        hg = HaploviewGenerator(
            encoded_result = encoded_result,
            locus_table = locus_table,
            ethnicity_info = eigen_result,
            output_path = "/the/output/folder")

    Generate the table::

        hg.haplotype_output(
            variant_list = HTN_variant_list,
            output_label = "HTN_1"
        )

    Where variant list should be in ``encoded_result.columns``.

    Two files will be saved under ``output_path``:

     - linkage_{output_label}.txt
     - map_{output_label}.map

    which can be used as the linkage format input for Haploview4.2
"""
import pandas as pd
from typing import List, Optional


class HaploviewGenerator():
    def __init__(
        self,
        encoded_result: pd.DataFrame,
        locus_table: pd.DataFrame,
        ethnicity_info: Optional[pd.DataFrame] = None,
        output_path: str = "/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/data_analysis_result/paper_1_output/linkage"
    ) -> None:
        self.encoded_result = encoded_result
        self.locus_table = locus_table
        if ethnicity_info is not None:
            self.ethnicity_info = ethnicity_info["ethnicity"]
        else:
            self.ethnicity_info = None
        self.output_path = output_path

    def encode(self, x: pd.Series) -> pd.DataFrame:
        genotype = x.name.split("-")[1]
        ref = genotype.split("/")[0]
        alt = genotype.split("/")[1]
        output = x.copy().replace(1, alt).replace(0, ref).replace(pd.NA, 0)
        output = pd.concat([output, output], axis = 1)
        output.columns = [f"{x.name}Allele1", f"{x.name}Allele2"]
        return output

    def haplotype_output(
        self,
        variant_list: List[str],
        output_label: str,
        ethnicity: Optional[str] = None
    ) -> None:
        """The API for generating output

        Args:
            variant_list: list[str], the list of variants to be used
            output_label: str, the name label for output under output_path
            ethnicity: Optional[str], the label used in ethnicity_info["ethnicity"],
                when specified, output will only include sample in this ethnicity
        """
        if (self.ethnicity_info is None) & (ethnicity is not None):
            raise ValueError("no ethnicity_info provided when initialize")
        elif (ethnicity is None):
            related_encoded_result = self.encoded_result
        else:
            related_encoded_result = self.encoded_result.loc[
                self.encoded_result.index.isin(
                    self.ethnicity_info.loc[self.ethnicity_info == ethnicity].index)]
        encoded_result_output = pd.DataFrame(data = related_encoded_result.index)
        encoded_result_output.columns = ["IndividualID"]
        encoded_result_output["PedID"] = 0
        encoded_result_output["FatherID"] = 0
        encoded_result_output["MotherID"] = 0
        encoded_result_output["Sex"] = 1
        encoded_result_output["Outcome"] = 0
        encoded_result_output = encoded_result_output[["PedID", "IndividualID", "FatherID", "MotherID", "Sex", "Outcome"]]
        for variant in variant_list:
            encoded_result_output = pd.merge(
                left = encoded_result_output,
                right = self.encode(related_encoded_result[variant]),
                how = "inner",
                left_on = "IndividualID",
                right_index = True)
        encoded_result_output.to_csv(
            f"{self.output_path}/linkage_{output_label}.txt",
            sep = "\t",
            header = False,
            index = False)
        self.locus_table.loc[
            [i.split("Allele")[0] for i in encoded_result_output.columns[6::2]],
            ["Pos"]
        ].to_csv(
            f"{self.output_path}/map_{output_label}.map",
            header = False,
            sep = "\t")
