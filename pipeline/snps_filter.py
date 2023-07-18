import pandas as pd

class SnpsFilter():
    """
    apply filter A and B to a one-hot SNPs table

    Example:
        filter_AB = SnpsFilter()
        filtered_result, drop_mask, drop_report = filter_AB.filter_A_and_B(df)

    Where df is a pd.DataFrame instance, with SNPs on the header, subject ID at the index,.
    i.e. This filter is dropping columns.
    """
    def __init__(self, drop_ratio: float = 0.1):
        self.drop_ratio = drop_ratio

    def filter_A(self, df, drop_ratio: float = None):
        """give filter_A's result: drop variants with >10% NA in population"""
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

        The value of filter be is based on unique number of non-NA values in each specfic position
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
