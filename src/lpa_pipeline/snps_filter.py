"""Give a first-step filter to all the SNPs encoded from Coassin's output

Apply the following:

    * filter A: drop SNPs with NA ratio higher than ``drop_ratio`` in population
    * filter B: drop SNPs with all 1's or all 0's (NA doesn't count for this)

Example::

    filter_AB = SnpsFilter(drop_ratio = 0.1)
    filtered_result, drop_mask, drop_report = filter_AB.filter_a_and_b(df)

Where df is a pd.DataFrame instance,
with SNPs ID on the header, subject ID at the index,
i.e. This filter is dropping columns.
"""
import pandas as pd
from typing import Optional, Tuple


class SnpsFilter:
    """
    apply filter A and B to a one-hot SNPs table

    Initialize:
        filter_AB = SnpsFilter(drop_ratio: float)

    Running:
        filtered_result, drop_mask, drop_report = filter_AB.filter_a_and_b(df)
    """

    def __init__(self, drop_ratio: float = 0.1) -> None:
        self.drop_ratio = drop_ratio

    def filter_a(
            self,
            df: pd.DataFrame,
            drop_ratio: Optional[float] = None
    ) -> pd.DataFrame:
        """give filter_A's result: drop variants with >drop_ratio NA in population"""
        if drop_ratio is None:
            drop_ratio = self.drop_ratio
        NA_count = df.isna().sum(axis=1)
        n_individual = df.shape[1]
        filtered_A = NA_count < n_individual * drop_ratio
        df = pd.concat([NA_count, filtered_A], axis=1)
        df.columns = ["NA_count", "filtered_A"]
        return df

    def filter_b(self, df: pd.DataFrame) -> pd.DataFrame:
        """add a filtered_B indicator on filter A result

        The value of filter B is based on unique number of non-NA values in each specific position
        """
        df = pd.DataFrame(df.nunique(axis=1) == 2, columns=["filtered_B"])
        return df

    def drop_report(
            self,
            df: pd.DataFrame,
            drop_ratio: Optional[float] = None
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Generate a drop report about filter A and B

        Args:
            df: pd.DataFrame, the table to be filtered
            drop_ratio: float, if provided, will cover the default in class initiation

        Returns:
            pd.DataFrame: the dropping boolean mask of df
            pd.DataFrame: the report for numbers of dropping
        """
        if drop_ratio is None:
            drop_ratio = self.drop_ratio
        df = df.T
        mask = pd.concat([self.filter_a(df),
                          self.filter_b(df)],
                         axis=1).drop(columns=["NA_count"])
        mask["filtered"] = mask["filtered_A"] & mask["filtered_B"]
        report = pd.DataFrame(mask[["filtered_A", "filtered_B", "filtered"]].sum(), columns=["left"])
        report["drop"] = df.shape[0] - report["left"]
        return mask, report

    def filter_a_and_b(
            self,
            df: pd.DataFrame,
            drop_ratio=None
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """ API run filter A and B and generate the result in one line

        Args:
            df: pd.DataFrame, the table to be filtered
            drop_ratio: float, if provided, will cover the default in class initiation

        Returns:
            pd.DataFrame: the table of SNPs passed the dropping
            pd.DataFrame: the dropping boolean mask of df
            pd.DataFrame: the report for numbers of dropping
        """
        drop_mask, drop_report = self.drop_report(df, drop_ratio=drop_ratio)
        filtered_df = df.T[drop_mask["filtered"]].T
        return filtered_df, drop_mask, drop_report


if __name__ == "__main__":
    import argparse
    import os
    import textwrap

    parser = argparse.ArgumentParser(prog="snps_filter.py",
                                     description=textwrap.dedent("""
        apply the following:
        filter A: drop SNPs with more than <drop_ratio> NA in population
        filter B: drop SNPs with all 1's or all 0's (NA doesn't count for this)

        Where the DataFrame passed in should be a csv file with SNPs on the header, subject ID at the index.
        i.e. This filter is dropping columns.
        """), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-I', "--input_path", type=str, required=True,
                        help="path to input encode result dataframe, in .csv format")
    parser.add_argument('-O', "--output_path", type=str, required=True,
                        help="path saving the output, this should be a folder rather than a file, "
                             "for 3 files will be saved")
    parser.add_argument('-r', '--drop_ratio', type=float, required=False, default=0.1,
                        help='the drop ratio for filter A, SNPs with NA ratio higher than this threshold will be '
                             'dropped, in float format, default 0.1 i.e. 10 percent')
    Args = parser.parse_args()

    read_in = pd.read_csv(Args.input_path, index_col=0)
    filter_AB = SnpsFilter(drop_ratio=Args.drop_ratio)
    filtered_result, filter_mask, filter_report = filter_AB.filter_a_and_b(read_in)
    filtered_result.to_csv(os.path.join(Args.output_path, "filtered_result.csv"))
    filter_mask.to_csv(os.path.join(Args.output_path, "drop_mask.csv"))
    filter_report.to_csv(os.path.join(Args.output_path, "drop_report.csv"))
