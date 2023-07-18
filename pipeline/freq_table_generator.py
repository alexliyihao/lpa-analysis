import pandas as pd

class FreqTableGenerator():
    """
    a generator class with API generate frequency table on all-one-hot variables,
    along with rarity check, originally designed for genome analysis

    Common usage:
    Given two pandas.DataFrame,
    class_info_table has one have columns <class_variable> that indicate each row belongs to,
    one_hot_table has all one-hot-variables,

    ftg = FreqTableGenerator()
    freq_table = ftg.generate_freq_table(
                    class_info_table = class_info_table,
                    one_hot_table = one_hot_table,
                    class_variable = "<class_variable>"
                    class_variable_list = ["<class_name_1>","<class_name_2>",...]
                    #if only need a part of <class_variable> column
                    )
    if need a rarity classification:
    freq_table_with_rarity = ftg.generate_freq_table_with_rarity(
                                class_info_table = class_info_table,
                                one_hot_table = one_hot_table,
                                class_variable = "<class_variable>"
                                class_variable_list = ["<class_name_1>","<class_name_2>",...]
                                #if only need a part of <class_variable> column
                                )
    """
    def __init__(self,
                 threshold:int = 0.01,
                 class_lower_case:bool = False,
                 encoding:dict = {0: "Not Detected",
                                  1: "Rare",
                                  2: "Common"}
                ):
        self.threshold = threshold
        self.class_lower_case = class_lower_case
        self.encoding = encoding

    def get_freq(self, class_string, class_wise_group,one_hot_table):
        class_subgroup = class_wise_group.get_group(class_string)
        class_subgroup, one_hot_table_subgroup = class_subgroup.align(one_hot_table, axis = 0, join = "inner")
        class_freq = pd.DataFrame({f"count_{class_string}":one_hot_table_subgroup.sum(axis = 0),
                                   f"total_{class_string}_detected":one_hot_table_subgroup.notna().sum(axis = 0),
                                   f"total_{class_string}_population":one_hot_table_subgroup.shape[0]
                                  })
        class_freq[f"freq_{class_string}"] = class_freq[f"count_{class_string}"] / class_freq[f"total_{class_string}_detected"]
        return class_freq

    def generate_freq_table(self, class_info_table,
                            one_hot_table,
                            class_variable: str,
                            class_variable_list = None,
                            return_class_variable_list = False):
        class_wise_group = class_info_table.groupby(class_variable)
        if class_variable_list is None:
            class_variable_list = class_wise_group.groups.keys()
        freq_list = [self.get_freq(
                        class_string = i,
                        class_wise_group = class_wise_group,
                        one_hot_table = one_hot_table)
                     for i in class_variable_list]
        freq_table = pd.concat(freq_list, axis = 1)
        if return_class_variable_list == True:
            return freq_table, class_variable_list
        else:
            return freq_table

    def common_check(self, freq_table,
                     class_variable_list:list,
                     threshold: float = None,
                     class_lower_case: bool = None,
                     encoding:dict = None):
        if threshold is None:
            threshold = self.threshold
        if class_lower_case is None:
            class_lower_case = self.class_lower_case
        if encoding is None:
            encoding = self.encoding

        for i in class_variable_list:
            if class_lower_case == True:
                i = i.lower()
            freq_name = f"freq_{i}"
            freq_table[i] = ((freq_table[freq_name] > 0) & (freq_table[freq_name] < 1)).astype(int)
            freq_table[i] = freq_table[i] + ((freq_table[freq_name] > threshold) &
                                             (freq_table[freq_name] < (1-threshold)))
            if encoding is not None:
                freq_table[i].replace(encoding, inplace = True)
        return freq_table

    def generate_freq_table_with_rarity(self,
                               class_info_table,
                               one_hot_table,
                               class_variable: str,
                               class_variable_list = None,
                               threshold: float = None,
                               class_lower_case: bool = None,
                               encoding: dict = None
                               ):
        freq_table, class_variable_list = self.generate_freq_table(
                                            class_info_table = class_info_table,
                                            one_hot_table = one_hot_table,
                                            class_variable = class_variable,
                                            class_variable_list = class_variable_list,
                                            return_class_variable_list = True)
        freq_table_with_rarity = self.common_check(
                                    freq_table = freq_table,
                                    class_variable_list = class_variable_list,
                                    threshold = threshold,
                                    class_lower_case = class_lower_case,
                                    encoding = encoding
                                    )
        return freq_table_with_rarity
