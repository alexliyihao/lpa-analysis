"""A generator computing relative frequency of SNP carrier by group

Common usage:

Given two pandas.DataFrame

* ``class_info_table`` has one have columns ``class_variable`` indicate the group
* ``one_hot_table`` has all one-hot-variables(SNPs encoding),

Compute the SNPs frequency in each class defined in <class_variable>:

    Initialize::

        ftg = freq_table_generator.FreqTableGenerator(
            threshold = 0.01
            encoding = {0: "Not Detected",
                        1: "Rare",
                        2: "Common"})

    Generate the table::

        freq_table = ftg.generate_freq_table(
            class_info_table = class_info_table,
            one_hot_table = one_hot_table,
            class_variable = "<class_variable>"
            class_variable_list = ["<class_name_1>","<class_name_2>",...]
            #if only need a part of <class_variable> column
            )

    If you need a rarity classification as columns as well::

        freq_table_with_rarity = ftg.generate_freq_table_with_rarity(
            class_info_table = class_info_table,
            one_hot_table = one_hot_table,
            class_variable = "<class_variable>"
            class_variable_list = ["<class_name_1>","<class_name_2>",...]
            #if only need a part of <class_variable> column
            )

Where threshold and encoding are for generating rarity:

* if the threshold < freq < 1- threshold, it will be encoded as 2
* if the 0 < freq < threshold or 1 - threshold < freq < 1, it will be encoded as 1
* if freq = 0 or freq = 1, it will be 0.

Thus, for encoding, please only modify the value, not the key.
"""
import pandas as pd


class FreqTableGenerator:
    """
    a generator class with API generate frequency table on all-one-hot variables,
    along with rarity check, originally designed for genome analysis

    Initialization:

        ftg = FreqTableGenerator(threshold: float,
            class_lower_case: bool,
            encoding: dict)

    Args:
        threshold: float, the threshold for rare vs. common for SNPs, default 0.01
        class_lower_case: bool, if True, all the class name will be in lower case
            for formatting
        encoding: dict{int: "str"}, the encoding for rarity, the keys should be
            kept as 0,1,2. 0 stands for there's no founding, 2 stands for the
            frequency is between threshold to 1-threshold
    """

    def __init__(self,
                 threshold: float = 0.01,
                 class_lower_case: bool = False,
                 encoding: dict = {0: "Not Detected",
                                   1: "Rare",
                                   2: "Common"}
                 ):
        self.threshold = threshold
        self.class_lower_case = class_lower_case
        self.encoding = encoding

    def get_freq(self, class_string: str, class_wise_group, one_hot_table):
        """Compute the freq for each class"""
        class_subgroup = class_wise_group.get_group(class_string)
        class_subgroup, one_hot_table_subgroup = class_subgroup.align(
            one_hot_table, axis=0, join="inner")
        class_freq = pd.DataFrame(
            {f"count_{class_string}": one_hot_table_subgroup.sum(axis=0),
             f"total_{class_string}_detected": one_hot_table_subgroup.notna().sum(axis=0),
             f"total_{class_string}_population": one_hot_table_subgroup.shape[0]})
        class_freq[f"freq_{class_string}"] = class_freq[f"count_{class_string}"] / class_freq[
            f"total_{class_string}_detected"]
        return class_freq

    def generate_freq_table(self,
                            class_info_table,
                            one_hot_table,
                            class_variable: str,
                            class_variable_list=None,
                            return_class_variable_list=False):
        """Generate the whole freq table for all class

        Args:
            class_info_table: pd.DataFrame, the table with class_variable
            one_hot_table: pd.DataFrame, the table with all the one-hot-encoded data
            class_variable: str, the name of column indicating the group setting
            class_variable_list: Optional[list[str]], if provided, only mentioned
                groups in class_variable will be listed
            return_class_variable_list: bool, if True, will return the
                class_variable_list use, even if not provided as variable

        Return:
            pd.DataFrame, the frequency table generated
            list[str], the list of group name included if return_class_variable_list == True
        """
        class_wise_group = class_info_table.groupby(class_variable)
        if class_variable_list is None:
            class_variable_list = class_wise_group.groups.keys()
        freq_list = [self.get_freq(
            class_string=i,
            class_wise_group=class_wise_group,
            one_hot_table=one_hot_table)
            for i in class_variable_list]
        freq_table = pd.concat(freq_list, axis=1)
        if return_class_variable_list:
            return freq_table, class_variable_list
        else:
            return freq_table

    def common_check(self, freq_table,
                     class_variable_list: list,
                     threshold: float = None,
                     class_lower_case: bool = None,
                     encoding: dict = None):
        """Append the commonness check to a frequency table

        Args:

            freq_table: pd.DataFrame, the frequency table provided
            class_variable_list: list[str], the list of group name included
            threshold: float, will cover the class setting if provided
            class_lower_case: bool, will cover the class setting if provided
            encoding: dict{int: "str"}, will cover the class setting if provided

        Returns:

            pd.DataFrame, the freq_table with rarity columns
        """
        if threshold is None:
            threshold = self.threshold
        if class_lower_case is None:
            class_lower_case = self.class_lower_case
        if encoding is None:
            encoding = self.encoding

        for i in class_variable_list:
            if class_lower_case:
                i = i.lower()
            freq_name = f"freq_{i}"
            freq_table[i] = ((freq_table[freq_name] > 0) & (freq_table[freq_name] < 1)).astype(int)
            freq_table[i] = freq_table[i] + ((freq_table[freq_name] > threshold) &
                                             (freq_table[freq_name] < (1 - threshold)))
            if encoding is not None:
                freq_table[i].replace(encoding, inplace=True)
        return freq_table

    def generate_freq_table_with_rarity(
            self,
            class_info_table,
            one_hot_table,
            class_variable: str,
            class_variable_list=None,
            threshold: float = None,
            class_lower_case: bool = None,
            encoding: dict = None
    ):
        """Generate the whole freq table for all class

        Args:
            class_info_table: pd.DataFrame, the table with class_variable
            one_hot_table: pd.DataFrame, the table with all the one-hot-encoded data
            class_variable: str, the name of column indicating the group setting
            class_variable_list: Optional[list[str]], if provided, only mentioned
                groups in class_variable will be listed
            threshold: float, will cover the class setting if provided
            class_lower_case: bool, will cover the class setting if provided
            encoding: dict{int: "str"}, will cover the class setting if provided

        Return:
            pd.DataFrame, the frequency table with rarity indicators
        """

        freq_table, class_variable_list = self.generate_freq_table(
            class_info_table=class_info_table,
            one_hot_table=one_hot_table,
            class_variable=class_variable,
            class_variable_list=class_variable_list,
            return_class_variable_list=True)

        freq_table_with_rarity = self.common_check(
            freq_table=freq_table,
            class_variable_list=class_variable_list,
            threshold=threshold,
            class_lower_case=class_lower_case,
            encoding=encoding)
        return freq_table_with_rarity


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="freq_table_generator.py",
                                     description="""
Given two table in csv format.
class_info_table has one have columns <class_variable> that indicate each row belongs to,
one_hot_table has all one-hot-encoded variables(SNPs)
Compute the SNPs' frequency in each class defined in <class_variable>

With a given <threshold> :
if the threshold < freq < 1- threshold, it will be encoded as "Common"
if the 0 < freq < threshold or 1 - threshold < freq < 1, it will be encoded as "Rare"
if freq = 0 or freq = 1, it will be "Not Detected".
""", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-C', "--class_info_table", type=str, required=True,
                        help="path to the table with class_info_data")
    parser.add_argument('-S', "--snps_table", type=str, required=True,
                        help="path to the table with one-hot-encode variable(SNPs)")
    parser.add_argument('-O', "--output", type=str, required=True,
                        help="path to the output, will be in csv format")
    parser.add_argument("-V", "--class_variable", type=str, required=True,
                        help="The class-indicating column in class_info_table")
    parser.add_argument("-L", "--class_variable_list", required=False, nargs='+',
                        help="If need a subset of class_variable, list as '-L EU AF HISP'")
    parser.add_argument('-T', "--threshold", type=float, required=False,
                        default=0.01, help="The threshold of Rare vs. Common, default 0.01")
    Args = parser.parse_args()

    ftg = FreqTableGenerator(threshold=Args.threshold)
    class_info_table = pd.read_csv(Args.class_info_table, index_col=0)
    snps_table = pd.read_csv(Args.snps_table, index_col=0)
    freq_table_with_rarity = ftg.generate_freq_table_with_rarity(
        class_info_table=class_info_table,
        one_hot_table=snps_table,
        class_variable=Args.class_variable,
        class_variable_list=Args.class_variable_list)
    freq_table_with_rarity.to_csv(Args.output)
