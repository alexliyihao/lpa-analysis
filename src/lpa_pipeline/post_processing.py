import pandas as pd
from typing import Optional
import statsmodels


class PostProcessor:
    """wrapper for all the postprocessing"""
    def __init__(
        self,
        locus_table: Optional[str] = None,
        p_value_threshold: float = 0.05,
        method: str = "fdr_bh"):
        self.cmd = CorrectMetalDirection(locus_table = locus_table)
        self.fa = FdrAdjustment(
            p_value_threshold = p_value_threshold,
            method = method)

    def post_process_association(self, df):
        df_output = df.reset_index()
        df_output = self.fa.appending_corrections_association(df = df_output)
        df_output = self.cmd.append_locus(df_output)
        return df_output.reset_index(drop = True)

    def post_process_meta_analysis(self, df):
        df_output = self.fa.appending_corrections_meta(df = df)
        df_output = self.cmd.correct_metal_complete(df = df_output)
        return df_output.reset_index(drop = True)

class CorrectMetalDirection:
    """correct the beta and direction discrepancy caused by metal"""

    def __init__(self, locus_table: str = None):
        self.locus_table = locus_table

    def correct_effect(self, x):
        """correct the effect direction, running for pd.DataFrame.apply()"""
        if x["metal_direction_correct"]:
            return(x["Effect"])
        else:
            return(x["Effect"] * (-1))

    def correct_direction(self, x):
        """correct the direction label, running for pd.DataFrame.apply()"""
        if x["metal_direction_correct"]:
            return(x["Direction"])
        else:
            return(x["Direction"].replace("+", "*").replace("-", "+").replace("*", "-"))

    def correct_metal_direction(self, df):
        """check each row for the metal direction and give corrected result"""
        assert((df["Variant"] == df["index"].str.extract("[0-9]+-[A|C|G|T]/([A/C/G/T])")[0]).all()),"Variant record incorrect"
        df["metal_direction_correct"] = (df["Variant"].str.lower() == df["Allele1"])
        df["Effect_corrected"] = df.apply(lambda x: self.correct_effect(x), axis = 1)
        df["Direction_corrected"] = df.apply(lambda x: self.correct_direction(x), axis = 1)
        return(df)

    def append_locus(self, df):
        """append locus table to the output"""
        if self.locus_table is None:
            raise ValueError("locus table is not given")
        df = df.loc[(df["index"] != "666-A/T")]
        df = pd.merge(
            left = df,
            right = self.locus_table,
            left_on = "index",
            right_index = True,
            how = "left"
            ).sort_values("FDR_adjusted_p-value")
        return df

    def correct_metal_complete(self, df):
        """complete pipeline correct the direction"""
        df = self.append_locus(df)
        df = self.correct_metal_direction(df)
        return df

class FdrAdjustment:
    "computing FDR adjustment for the output of association or metal_toolkit"
    def __init__(
        self,
        p_value_threshold: float = 0.05,
        method: str = "fdr_bh"):
        self.p_value_threshold = p_value_threshold
        self.method = method

    def appending_corrections_meta(
        self,
        df: pd.DataFrame
        ) -> pd.DataFrame:
        """compute the FDR corrected p-value for each trait in meta-analysis result

        Arg:
            df: pd.DataFrame, the output table from meta-analysis
                               module
        Return:
            pd.DataFrame, df with FDR, Sidak and Bonferonni info
        """
        output = []
        for trait in df["trait"].unique():
            df_sub = df.loc[
                (df["trait"] == trait) &
                (df["P-value"].notna())]

            reject, pvals_corrected, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(
                pvals = df_sub["P-value"],
                alpha = self.p_value_threshold,
                method = self.method
            )
            df_sub["FDR_adjusted_p-value"] = pvals_corrected
            df_sub["significant_FDR"] = reject
            df_sub["Sidak_adjusted_alpha"] = alphacSidak
            df_sub["Bonferonni_adjusted_alpha"] = alphacBonf
            output.append(df_sub)
        return pd.concat(output, axis = 0).sort_values("FDR_adjusted_p-value")

    def appending_corrections_association(
        self,
        df: pd.DataFrame) -> pd.DataFrame:
        """compute the FDR corrected p-value for assocation result

        Arg:
            df: pd.DataFrame, the output table from association module
        Return:
            pd.DataFrame, df with FDR, Sidak and Bonferonni info
        """
        reject, pvals_corrected, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(
            pvals = df["snp_pos P-value"],
            alpha = self.p_value_threshold,
            method = self.method
        )
        df["FDR_adjusted_p-value"] = pvals_corrected
        df["significant_FDR"] = reject
        df["Sidak_adjusted_alpha"] = alphacSidak
        df["Bonferonni_adjusted_alpha"] = alphacBonf
        return df.sort_values("FDR_adjusted_p-value")

    def appending_meta_analyis(self, path: str):
        """wrapper for appending_corrections_meta"""
        df = pd.read_excel(path)
        self.appending_corrections_meta(df).to_excel(
            path.replace(".xlsx","_with_corrected_significance.xlsx")
            )

    def appending_association(self, path: str):
        """appending_corrections_association"""
        df = pd.read_csv(path)
        self.appending_corrections_association(df).reset_index(
            drop = True
        ).to_excel(path.replace(".csv","_with_corrected_significance.xlsx"))


class CorrectMetalDirection:
    """correct the beta and direction discrepancy caused by metal"""

    def __init__(self, locus_table: str = None):
        self.locus_table = locus_table

    def correct_effect(self, x):
        """correct the effect direction, running for pd.DataFrame.apply()"""
        if x["metal_direction_correct"]:
            return(x["Effect"])
        else:
            return(x["Effect"] * (-1))

    def correct_direction(self, x):
        """correct the direction label, running for pd.DataFrame.apply()"""
        if x["metal_direction_correct"]:
            return(x["Direction"])
        else:
            return(x["Direction"].replace("+", "*").replace("-", "+").replace("*", "-"))

    def correct_metal_direction(self, df):
        """check each row for the metal direction and give corrected result"""
        assert((df["Variant"] == df["index"].str.extract("[0-9]+-[A|C|G|T]/([A/C/G/T])")[0]).all()),"Variant record incorrect"
        df["metal_direction_correct"] = (df["Variant"].str.lower() == df["Allele1"])
        df["Effect_corrected"] = df.apply(lambda x: self.correct_effect(x), axis = 1)
        df["Direction_corrected"] = df.apply(lambda x: self.correct_direction(x), axis = 1)
        return(df)

    def append_locus(self, df):
        """append locus table to the output"""
        if self.locus_table is None:
            raise ValueError("locus table is not given")
        df = pd.merge(
            left = df.loc[(df["index"] != "666-A/T")],
            right = self.locus_table,
            left_on = "index",
            right_index = True,
            how = "left"
            ).sort_values("FDR_adjusted_p-value")
        return df

    def correct_metal_complete(self, df):
        """complete pipeline correct the direction"""
        df = self.append_locus(df)
        df = self.correct_metal_direction(df)
        return df

class FdrAdjustment:
    "computing FDR adjustment for the output of association or metal_toolkit"
    def __init__(
        self,
        p_value_threshold: float = 0.05,
        method: str = "fdr_bh"):
        self.p_value_threshold = p_value_threshold
        self.method = method

    def appending_corrections_meta(
        self,
        df: pd.DataFrame
        ) -> pd.DataFrame:
        """compute the FDR corrected p-value for each trait in meta-analysis result

        Arg:
            df: pd.DataFrame, the output table from meta-analysis
                               module
        Return:
            pd.DataFrame, df with FDR, Sidak and Bonferonni info
        """
        output = []
        for trait in df["trait"].unique():
            df_sub = df.loc[
                (df["trait"] == trait) &
                (df["P-value"].notna())]

            reject, pvals_corrected, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(
                pvals = df_sub["P-value"],
                alpha = self.p_value_threshold,
                method = self.method
            )
            df_sub["FDR_adjusted_p-value"] = pvals_corrected
            df_sub["significant_FDR"] = reject
            df_sub["Sidak_adjusted_alpha"] = alphacSidak
            df_sub["Bonferonni_adjusted_alpha"] = alphacBonf
            output.append(df_sub)
        return pd.concat(output, axis = 0).sort_values("FDR_adjusted_p-value")

    def appending_corrections_association(
        self,
        df: pd.DataFrame) -> pd.DataFrame:
        """compute the FDR corrected p-value for assocation result

        Arg:
            df: pd.DataFrame, the output table from association module
        Return:
            pd.DataFrame, df with FDR, Sidak and Bonferonni info
        """
        reject, pvals_corrected, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(
            pvals = df["snp_pos P-value"],
            alpha = self.p_value_threshold,
            method = self.method
        )
        df["FDR_adjusted_p-value"] = pvals_corrected
        df["significant_FDR"] = reject
        df["Sidak_adjusted_alpha"] = alphacSidak
        df["Bonferonni_adjusted_alpha"] = alphacBonf
        return df.sort_values("FDR_adjusted_p-value")

    def appending_meta_analyis(self, path: str):
        """wrapper for appending_corrections_meta"""
        df = pd.read_excel(path)
        self.appending_corrections_meta(df).to_excel(
            path.replace(".xlsx","_with_corrected_significance.xlsx")
            )

    def appending_association(self, path: str):
        """appending_corrections_association"""
        df = pd.read_csv(path)
        self.appending_corrections_association(df).reset_index(
            drop = True
        ).to_excel(path.replace(".csv","_with_corrected_significance.xlsx"))
