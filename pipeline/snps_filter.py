"""
apply the following:
filter A: drop SNPs with more than <drop_ratio> NA in population
filter B: drop SNPs with all 1's or all 0's (NA doesn't count for this)

Example:
    filter_AB = SnpsFilter(drop_ratio = 0.1)
    filtered_result, drop_mask, drop_report = filter_AB.filter_A_and_B(df)

Where df is a pd.DataFrame instance, with SNPs on the header, subject ID at the index,.
i.e. This filter is dropping columns.
"""
import pandas as pd

class SnpsFilter():
    """
    apply filter A and B to a one-hot SNPs table

    Initialize:
        filter_AB = SnpsFilter(drop_ratio: float)

    Running:
        filtered_result, drop_mask, drop_report = filter_AB.filter_A_and_B(df)
    """
    def __init__(self, drop_ratio: float = 0.1):
        self.drop_ratio = drop_ratio

    def filter_A(self, df, drop_ratio: float = None):
        """give filter_A's result: drop variants with >drop_ratio NA in population"""
        if drop_ratio is None:
            drop_ratio = self.drop_ratio
        NA_count = df.isna().sum(axis = 1)
        n_individual = df.shape[1]
        filtered_A =  NA_count < n_individual * drop_ratio
        df = pd.concat([NA_count,filtered_A], axis =1)
        df.columns = ["NA_count", "filtered_A"]
        return df

    def filter_B(self, df):
        """add a filtered_B indicator on filter A result

        The value of filter B is based on unique number of non-NA values in each specfic position
        """
        df = pd.DataFrame(df.nunique(axis = 1) == 2, columns = ["filtered_B"])
        return df

    def drop_report(self, df, drop_ratio = None):
        if drop_ratio is None:
            drop_ratio = self.drop_ratio
        df = df.T
        mask = pd.concat([self.filter_A(df),
                          self.filter_B(df)],
                        axis = 1).drop(columns=["NA_count"])
        mask["filtered"] = mask["filtered_A"] & mask["filtered_B"]
        report = pd.DataFrame(mask[["filtered_A","filtered_B", "filtered"]].sum(), columns = ["left"])
        report["drop"] = df.shape[0] - report["left"]
        return (mask,report)

    def filter_A_and_B(self, df, drop_ratio = None):
        drop_mask, drop_report = self.drop_report(df, drop_ratio = drop_ratio)
        filtered_df = df.T[drop_mask["filtered"]].T
        return(filtered_df, drop_mask, drop_report)

if __name__ == "__main__":
    import argparse
    import os
    parser = argparse.ArgumentParser(prog = "snps_filter.py",
                                     description=\
"""
apply the following:
filter A: drop SNPs with more than <drop_ratio> NA in population
filter B: drop SNPs with all 1's or all 0's (NA doesn't count for this)

Where the DataFrame passed in should be a csv file with SNPs on the header, subject ID at the index.
i.e. This filter is dropping columns.
""",formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-I', "--input_path", type = str, required = True,
                        help = "path to input encode result dataframe, in .csv format")
    parser.add_argument('-O', "--output_path", type = str, required = True,
                        help = "path saving the output, this should be a folder rather than a file, for 3 files will be saved")
    parser.add_argument('-r', '--drop_ratio', type = float, required = False, default = 0.1,
                        help = 'the drop ratio for filter A, SNPs with NA ratio higher than this threshold will be dropped, in float format, default 0.1 i.e. 10 percent')
    Args = parser.parse_args()

    df = pd.read_csv(Args.input_path, index_col = 0)
    filter_AB = SnpsFilter(drop_ratio = Args.drop_ratio)
    filtered_result, drop_mask, drop_report = filter_AB.filter_A_and_B(df)
    filtered_result.to_csv(os.path.join(Args.output_path, "filtered_result.csv"))
    drop_mask.to_csv(os.path.join(Args.output_path, "drop_mask.csv"))
    drop_report.to_csv(os.path.join(Args.output_path, "drop_report.csv"))
