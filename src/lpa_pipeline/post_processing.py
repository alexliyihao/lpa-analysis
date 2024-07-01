"""Append necessary data/correction for association/meta-analysis outputs

For each association/meta-analysis output, the module does the following:

* Compute the FDR-adjusted p-value from ``statsmodels.stats.multitest.multipletests``,
* Correct the ``METAL``'s direction discrepancy,
  ``METAL`` is not computing for the direction of alternative but the first allele it met
* the ``post_processing`` module also appended the locus information necessary for analysis

Example::

    pp = PostProcessor(locus_table = locus_table,snp_alias = "variant")
    pp.post_process_association(df = some_association_result)
    pp.post_process_meta_analysis(df = aggregate_result)

"""
import pandas as pd
from typing import Optional
import statsmodels
import statsmodels.stats.multitest


class CorrectMetalDirection:
    """correct the beta and direction discrepancy caused by metal"""

    def __init__(self, locus_table: pd.DataFrame = None) -> None:
        self.locus_table = locus_table

    @staticmethod
    def correct_effect(x: pd.DataFrame) -> pd.DataFrame:
        """correct the effect direction, running for pd.DataFrame.apply()"""
        if x["metal_direction_correct"]:
            return x["Effect"]
        else:
            return x["Effect"] * (-1)

    @staticmethod
    def correct_direction(x: pd.DataFrame) -> pd.DataFrame:
        """correct the direction label, running for pd.DataFrame.apply()"""
        if x["metal_direction_correct"]:
            return x["Direction"]
        else:
            return x["Direction"].replace("+", "*").replace("-", "+").replace("*", "-")

    def correct_metal_direction(self, df: pd.DataFrame, mode: str = "variant") -> pd.DataFrame:
        """check each row for the metal direction and give corrected result"""
        if mode == "variant":
            assert ((df["Alt"] == df["index"].str.extract("[0-9]+-[A|C|G|T]/([A/C/G/T])")[
                0]).all()), "Variant record incorrect"
        df["metal_direction_correct"] = (df["Alt"].str.lower() == df["Allele1"])
        df["Effect_corrected"] = df.apply(lambda x: self.correct_effect(x), axis=1)
        df["Direction_corrected"] = df.apply(lambda x: self.correct_direction(x), axis=1)
        return df

    def append_locus(self, df: pd.DataFrame) -> pd.DataFrame:
        """append locus table to the output

        Args:
            df: pd.DataFrame, the result from either association or meta-analysis

        Returns:
            pd.DataFrame: df merged with the locus table
        """
        if self.locus_table is None:
            raise ValueError("locus table is not given")
        #df = df.loc[(df["index"] != "666-A/T") & (df["index"] != "594-A/G") & (df["index"] != "621-T/C")]
        df = pd.merge(
            left=df,
            right=self.locus_table,
            left_on="index",
            right_index=True,
            how="left"
        ).sort_values("FDR_adjusted_p-value")
        return df

    def correct_metal_complete(self, df: pd.DataFrame, mode: Optional[str] = "variant") -> pd.DataFrame:
        """complete pipeline correct the direction

        Args:
            df: pd.DataFrame, the result from meta-analysis (metal_toolkit module)
            mode: Optional[str], if using "variant", will check the consistency between results
        Result:
            pd.DataFrame, the result with corrected effect and corrected direction column
        """
        df = self.append_locus(df)
        df = self.correct_metal_direction(df, mode=mode)
        return df


class FdrAdjustment:
    """computing FDR adjustment for the output of association or metal_toolkit"""

    def __init__(
            self,
            p_value_threshold: float = 0.05,
            method: str = "fdr_bh",
            snp_alias: str = "variant"
    ) -> None:
        self.p_value_threshold = p_value_threshold
        self.method = method
        self.snp_alias = snp_alias

    def appending_corrections_meta(
            self,
            df: pd.DataFrame
    ) -> pd.DataFrame:
        """compute the FDR corrected p-value for each trait in meta-analysis result

        Args:
            df: pd.DataFrame, the output table from meta-analysis
                               module
        Returns:
            pd.DataFrame, df with FDR, Sidak and Bonferonni info
        """
        output = []
        for trait in df["trait"].unique():
            df_sub = df.loc[
                (df["trait"] == trait) &
                (df["P-value"].notna())]

            reject, pvals_corrected, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(
                pvals=df_sub["P-value"],
                alpha=self.p_value_threshold,
                method=self.method
            )
            df_sub["FDR_adjusted_p-value"] = pvals_corrected
            df_sub["significant_FDR"] = reject
            df_sub["Sidak_adjusted_alpha"] = alphacSidak
            df_sub["Bonferonni_adjusted_alpha"] = alphacBonf
            output.append(df_sub)
        return pd.concat(output, axis=0).sort_values("FDR_adjusted_p-value")

    def appending_corrections_association(
            self,
            df: pd.DataFrame
    ) -> pd.DataFrame:
        """compute the FDR corrected p-value for association result

        Args:
            df: pd.DataFrame, the output table from association module
        Returns:
            pd.DataFrame, df with FDR, Sidak and Bonferonni info
        """
        reject, pvals_corrected, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(
            pvals=df[f"{self.snp_alias} P-value"],
            alpha=self.p_value_threshold,
            method=self.method
        )
        df["FDR_adjusted_p-value"] = pvals_corrected
        df["significant_FDR"] = reject
        df["Sidak_adjusted_alpha"] = alphacSidak
        df["Bonferonni_adjusted_alpha"] = alphacBonf
        return df.sort_values("FDR_adjusted_p-value")

    def appending_meta_analysis(self, path: str) -> None:
        """wrapper for appending_corrections_meta"""
        df = pd.read_excel(path)
        self.appending_corrections_meta(df).to_excel(
            path.replace(".xlsx", "_with_corrected_significance.xlsx")
        )

    def appending_association(self, path: str) -> None:
        """appending_corrections_association"""
        df = pd.read_csv(path)
        self.appending_corrections_association(df).reset_index(
            drop=True
        ).to_excel(path.replace(".csv", "_with_corrected_significance.xlsx"))


class PostProcessor:
    """wrapper for all the postprocessing"""

    def __init__(
            self,
            locus_table: Optional[str] = None,
            p_value_threshold: Optional[float] = 0.05,
            method: Optional[str] = "fdr_bh",
            snp_alias: Optional[str] = "variant",
            ethnicity: Optional[list] = ["EU", "AF", "HISP"]
    ) -> None:
    """initializer of post-processing procedure

    Args:
        locus_table: Optional[pd.DataFrame], the locus table generate by locus_collector module
        p_value_threshold: Optional[float], used for the FDR-adjusted p-value, default 0.05
        method: Optional[str], the method argument used for statsmodels.stats.multitest.multipletests
        snp_alias: Optional[str], the alias of the variable of interest, usually unified with the
            snp_alias setting in association and metal_toolkit modules,
        ethnicity: Optional[list], the ethnicity label
    """
        self.cmd = CorrectMetalDirection(locus_table=locus_table)
        self.fa = FdrAdjustment(
            p_value_threshold=p_value_threshold,
            method=method,
            snp_alias=snp_alias)
        self.ethnicity = ethnicity
        self.snp_alias = snp_alias
        self.tidy_output = [
            'index',
            'trait', 'Effect_corrected', 'Direction_corrected',
            'FDR_adjusted_p-value', 'significant_FDR',
            'mylocus', 'coding', 'wt', 'mut', 'novel',
            'total_count', 'total_population']

    def post_process_association(self, df: pd.DataFrame) -> pd.DataFrame:
        """The wrapper of all post process steps for an association module result

        It computes FDR-adjusted p-value and append locus information for the result

        Args:
            df: pd.DataFrame, the output of association module

        Returns:
            pd.DataFrame: the module with FDR adjusted p-value and locus info
        """
        df_output = df.reset_index()
        df_output = self.fa.appending_corrections_association(df=df_output)
        df_output = self.cmd.append_locus(df_output)
        return df_output.reset_index(drop=True)

    def post_process_meta_analysis(self, df: pd.DataFrame, mode: str = "variant") -> pd.DataFrame:
        """The wrapper of all post process steps for a metal_toolkit module result

        It computes FDR-adjusted p-value, correct the direction discrepancy from METAL,
        and append locus information for the result

        Args:
            df: pd.DataFrame, the output of association module

        Returns:
            pd.DataFrame: the module with FDR adjusted p-value, correct direction and locus info
        """
        df_output = self.fa.appending_corrections_meta(df=df)
        df_output = self.cmd.correct_metal_complete(df=df_output, mode=mode)
        return df_output.reset_index(drop=True)

    def clean_output(self, df: pd.DataFrame) -> pd.DataFrame:
        """An API to reorder the output to better readability

        Args:
            df: pd.DataFrame, the output need to be cleaned

        Returns:
            pd.DataFrame: the format-cleaned output
        """
        ethnicity_columns: list = [
            [f'{self.snp_alias} Beta_{i}',
             f'{self.snp_alias} P-value_{i}',
             f'rel_freqs_{i}',
             f'abs_freqs_{i}',
             f'n_sample_{i}']
            for i in self.ethnicity]
        tidied_columns: list = self.tidy_output + [
                column for columns in ethnicity_columns for column in columns]
        return df.loc[:, tidied_columns]
