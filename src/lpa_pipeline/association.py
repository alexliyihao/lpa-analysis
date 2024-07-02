"""This pipeline running an iterative association test for each snps.

The iterative strategy is as follows:

* for each the target variable in ``target_strategy``
    #. Extract the column from ``target_dataset``:
    #. run preprocessing, if the value of key ``preprocessing`` is given in ``target_strategy``
    #. group the dataset by columns defined by ``extra_iterate_on`` if provided
    #. for each group, iterate over columns of ``encoded_snp`` (each SNP):
        1. concatenate it with ``other_exogs``, generate the exogenous dataset
        2. if ``one_to_one_exogs`` is provided, use ``one_to_one_strategy``
           finding other columns and concatenate them as well
        3. run regressions specified by ``target_strategy``.engine
        4. combine the results from each SNP
        5. save the regression output

Two APIs provided
 * sklearn style:

     * 3-line style::

        snp_asso = SNPAssociation()
        snp_asso.fit(**kwargs)
        snp_asso.transform()

     * 2-line style::

        snp_asso = SNPAssociation()
        snp_asso.fit_transform(**kwargs)

 * function style::

        snp_asso = SNPAssociation()
        snp_asso.association_test(**kwargs_1) #kwargs_1 can be new kwargs
"""

import gc
import json
import os
import warnings
from typing import Dict, List, Optional, Union, Callable, Tuple, Hashable
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from statsmodels.discrete.discrete_model import BinaryResultsWrapper
from statsmodels.regression.linear_model import RegressionResultsWrapper
from tqdm import tqdm

# ----------------------------------Encoding Strategy Example--------------------------------------

def encode_dementia(df: pd.DataFrame) -> pd.DataFrame:
    """An example preprocessing strategy used in target_strategy,

    It should be a callable takes in a pd.DataFrame and returns another one
    """
    df.dropna(inplace=True)
    df = df[df != 3] - 1
    return df


def dropna(df: pd.DataFrame) -> pd.DataFrame:
    """A simpler preprocessing strategy used in target_strategy

    It should be a callable takes in a pd.DataFrame and returns another one
    """
    return df.dropna()


def filter_c(df: pd.DataFrame, threshold=5) -> pd.DataFrame:
    """The filter C on existing frequency,

    both 0 and 1 should appear for more than <threshold> times

    Args:
        df: pandas DataFrame, the index should be the Sample ID, and the columns are the snps
        threshold: int, the number threshold that both 0 and 1 should appear more than this number
    """
    filtered_snp = [i for i in df.columns if sum(df[i].value_counts() > threshold) > 1]
    return df[filtered_snp]


def target_strategy() -> Dict[str, Dict[str, Union[sm.Logit, sm.OLS, Callable]]]:
    """
    An example target_strategy used in the association analysis

    The key for main dict is the variable name
    The value is a dict with the following keys:
      - "engine", whose value is a statsmodels model
      - "preprocessing", whose value is a callable takes in a pd.DataFrame and returns another one
    """
    return {'STROKE': {"engine": sm.Logit,
                       "preprocessing": dropna},
            'DEMENTIA': {"engine": sm.Logit,
                         "preprocessing": encode_dementia},
            'DIABETES': {"engine": sm.Logit,
                         "preprocessing": dropna},
            'HEART': {"engine": sm.Logit,
                      "preprocessing": dropna},
            "HYPERTENSION": {"engine": sm.Logit,
                             "preprocessing": dropna},
            'LIP01_B': {"engine": sm.OLS,
                        "preprocessing": dropna},
            'LIP02_B': {"engine": sm.OLS,
                        "preprocessing": dropna},
            'LIP03_B': {"engine": sm.OLS,
                        "preprocessing": dropna},
            'LIP04_2': {"engine": sm.OLS,
                        "preprocessing": dropna},
            'INSL01_2': {"engine": sm.OLS,
                         "preprocessing": dropna},
            'INSL02_2': {"engine": sm.OLS,
                         "preprocessing": dropna},
            'INSL03_2': {"engine": sm.OLS,
                         "preprocessing": dropna},
            'HBA1C_2': {"engine": sm.OLS,
                        "preprocessing": dropna}
            }


def target_strategy_serum() -> Dict[str, Dict[str, Union[sm.Logit, sm.OLS, Callable]]]:
    """
    An example target_strategy used in the serum analysis
    """
    return {'lpa': {"engine": sm.OLS,
                    "preprocessing": None},
            'wIS': {"engine": sm.OLS,
                    "preprocessing": None}
            }


# ----------------------------------Iterator API--------------------------------------

class SNPAssociation:
    """a class for running SNP association pipeline"""

    def __init__(self) -> None:
        pd.set_option('mode.chained_assignment', None)
        warnings.simplefilter('ignore', ConvergenceWarning)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        self._encoded_snp = None
        self._other_exogs = None
        self._target_dataset = None
        self._target_strategy = None
        self._output_path = None
        self._snp_alias = None
        self._one_to_one_exogs = None
        self._one_to_one_strategy = None
        self._one_to_one_alias = None
        self._na_strategy = None
        self._group_na_strategy = None
        self._extra_iterate_on = None
        self._snps_preprocessing_strategy = None
        self._verbose = None

    def debug_log(self, *args) -> None:
        """only print the result when it's in verbose mode == 2"""
        if self._verbose == 2:
            print(*args)

    def fit(self,
            encoded_snp: pd.DataFrame,
            other_exogs: pd.DataFrame,
            target_dataset: pd.DataFrame,
            target_strategy: Dict[str, Dict[str, Union[sm.Logit, sm.OLS, Callable]]],
            output_path: str,
            snp_alias: str = "variant",
            one_to_one_exogs: Optional[pd.DataFrame] = None,
            one_to_one_strategy: Optional[Callable] = None,
            one_to_one_alias: Optional[str] = None,
            na_strategy: str = "drop",
            group_na_strategy: str = "snp_wise",
            extra_iterate_on: Optional[list] = None,
            snps_preprocessing_strategy: Optional[Callable] = None,
            verbose: int = 0
            ) -> None:
        """API setting up the regression (not running)

        Args:

            encoded_snp: pd.DataFrame, the dataframe to be looped on columns

            other_exogs: pd.DataFrame, the dataframe taking all the other variables

            target_dataset: pd.DataFrame, the dataframe taking all the target variables

            target_strategy: dict[str, dict[str, funcs or models]],
                 The dictionary provide pre-processing to specific column,
                 Only column mentioned in keys will be included in the running.
                 The inner dictionary should have two keys:

                  - "engine": statsmodels.api models
                     designed for statsmodels.discrete.discrete_model.Logit or
                     statsmodels.regression.linear_model.OLS,
                     but any model's .fit results provide .params .bse .pvalues will work

                  - "preprocessing": funcs
                     the function should take a pd.DataFrame/pd.Series as input
                     and a pd.DataFrame/pd.Series as output,
                     If this strategy is None, the column will be used as
                     exogenous variables as-is

            output_path: str, the output root path,
                         the saving name will be printed for reference

            snp_alias: Optional[str], the name used for snp column in output DataFrame

            one_to_one_exogs: Optional[pd.DataFrame], the dataframe providing
                              specific exog variable based on encoded_snp

            one_to_one_strategy: Optional[funcs], when given a snp name, this
                                 function should output the corresponding
                                 one_to_one_exogs column name

            one_to_one_alias: Optional[str]: the name used for one-to-one variables
                              in output DataFrame

            na_strategy: Optional[str], the strategy statsmodels dealing with NAs,
                         Available options are ‘none’, ‘drop’, and ‘raise’.
                         If ‘none’, no nan checking is done.
                         If ‘drop’, any observations with nans are dropped.
                         If ‘raise’, an error is raised.

            group_na_strategy: Optional[str], how na_strategy apply when have
                               extra_iterate_on, can be "snp-wise" or "group-wise",
                               it's working for status when na_strategy == "drop"

            extra_iterate_on: Optional[list], the list of extra iterate variable in other_exog
                              please be noticed that it may cause severe explosions in running time

            snps_preprocessing_strategy: Optional[funcs], the preprocessing applied to snps
                           before the regression analysis
                           this function should take a pd.DataFrame as input and output

            verbose: Optional[int] in {0,1,2} default 0
                     if verbose = 1, the regression will give the saving path
                     if verbose = 2, the regression will output all the related
                     values during the computation, only for debugging purpose,
                     use with care for it will create massive I/O
        """
        self._encoded_snp = encoded_snp
        self._other_exogs = other_exogs
        self._target_dataset = target_dataset
        self._target_strategy = target_strategy
        self._output_path = output_path
        self._snp_alias = snp_alias
        self._one_to_one_exogs = one_to_one_exogs
        self._one_to_one_strategy = one_to_one_strategy
        self._one_to_one_alias = one_to_one_alias
        self._na_strategy = na_strategy
        self._group_na_strategy = group_na_strategy
        self._extra_iterate_on = extra_iterate_on
        self._snps_preprocessing_strategy = snps_preprocessing_strategy
        self._verbose = verbose

    def transform(self, output_path: str = None, verbose: int = None) -> None:
        """API Run the actual association test

        Args:
            output_path: Optional[str]: if provided, only in this run, it will cover
                         the output_path previously provided
            verbose: Optional[int]: if provided, only in this run, it will cover
                         the verbose previously provided
        """
        if output_path is None:
            output_path = self._output_path
        if verbose is None:
            verbose = self._verbose

        return self.association_test(
            encoded_snp=self._encoded_snp,
            other_exogs=self._other_exogs,
            target_dataset=self._target_dataset,
            target_strategy=self._target_strategy,
            output_path=output_path,
            one_to_one_exogs=self._one_to_one_exogs,
            one_to_one_strategy=self._one_to_one_strategy,
            one_to_one_alias=self._one_to_one_alias,
            na_strategy=self._na_strategy,
            extra_iterate_on=self._extra_iterate_on,
            group_na_strategy=self._group_na_strategy,
            snps_preprocessing_strategy=self._snps_preprocessing_strategy,
            snp_alias=self._snp_alias,
            verbose=verbose
        )

    def fit_transform(self, **kwargs):
        """API run fit and transform in one method"""
        self.fit(**kwargs)
        self.transform()

    def _get_column_name(self, col: List[str]) -> List[str]:
        """for a list of given column, generate suffix necessary as header"""
        col_set = [[f"{i} Beta", f"{i} SE", f"{i} P-value"] for i in col]
        return [i for j in col_set for i in j]

    def _generate_filename(
            self,
            n_snp: int,
            endog_variable: str,
            engine: str,
            extra_label: Optional[Tuple[Tuple[str], Tuple[str]]]) -> str:
        """generate a file name"""
        if extra_label is None:
            extra = ""
        else:
            extra = "_".join([f"{case[0]}={case[1]}"
                              for case in zip(extra_label[0], extra_label[1])]) + "_"
        return f"{endog_variable}_{engine}_{extra}N_snp={n_snp}.csv"

    def _generate_info(
            self,
            endog: Union[pd.Series, pd.DataFrame],
            engine: str,
            variable_name: List[str],
            extra_label: Optional[Tuple[Tuple[str], Tuple[str]]],
            processed_snp: List[str],
            n_sample: int,
            errors: Dict[str, str]) -> Dict[str, Union[str, int, List[str]]]:
        """generate an info dictionary for a regression task"""
        if isinstance(endog, pd.core.series.Series):
            endog_variable = endog.name
        elif isinstance(endog, pd.core.frame.DataFrame):
            endog_variable = endog.columns[0]
        else:
            raise ValueError(
                "not a proper endogenous variable format, "
                "it should be either a pd.Series or pd.DataFrame"
            )
        N_snp = len(processed_snp)
        return {"filename": self._generate_filename(
            n_snp=N_snp,
            endog_variable=endog_variable,
            engine=engine,
            extra_label=extra_label),
            "formula": f"{endog_variable} ~ {' + '.join(variable_name)}",
            "N_sample": n_sample,
            "N_snp": N_snp,
            "errors": errors
        }

    def _save_association_result(
            self,
            result: pd.DataFrame,
            info: Dict[str, Union[str, int, List[str]]],
            output_path: str,
            verbose: int):
        """save result and info to output path as the name specified in info"""
        # update the filename
        saving_path = os.path.join(output_path, info["filename"])
        info["filename"] = saving_path
        # save the csv file
        result.to_csv(saving_path)
        # save a txt file with the same name
        with open(saving_path.replace("csv", "txt"), "w+") as file:
            # write info into this file
            json.dump(info, file, indent=4)
        if verbose == 1:
            print(f"{saving_path}/.txt saved")

    def _calculate_freq(self, df: pd.DataFrame, cols: str) -> Tuple[float, int, int]:
        """
        given a pandas.dataFrame <df> and specific column name <cols>,
        compute the frequency

        Please note that this is only designed for 0/1 or boolean column
        When using quantitative value the freq it's not making realistic sense,
        So the computing is awkward... just want to keep the code running even
        with quantitative columns. Feel free to override it.

        Args:
            df: pd.DataFrame, a data table
            cols: str, column name in df

        Returns:
            tuple containing

            - float: the relative frequency of 1s in this column.
            - float: the number of 1s in the column.
            - float: the number of rows in the column.
        """
        col: pd.Series = df[cols]
        return (col.sum() / len(col)), col.sum(), len(col)

    def _association_snp(
            self,
            snp_table: pd.DataFrame,
            other_exogs: pd.DataFrame,
            endog: Union[pd.DataFrame, pd.Series],
            engine: Callable,
            one_to_one_exogs: Optional[pd.DataFrame] = None,
            one_to_one_strategy: Optional[Callable] = None,
            one_to_one_alias: List[str] = None,
            snps_preprocessing_strategy: Optional[Callable] = None,
            na_strategy: str = "drop",
            group_na_strategy: str = "snp_wise",
            snp_alias: str = "variant",
            extra_label: Optional[Tuple[List[str], List[Hashable]]] = None
    ) -> Tuple[pd.DataFrame, Dict]:
        """association test designed for snp variants analysis

        For each column in snp_table, run an individual <endog>~column+<other> on <engine>
        return the beta/weight/effect size, standard error and p-values

        Args:

            snp_table: pd.DataFrame, the table encoding snps,
                       the columns is snp name, rows are individuals, values are numerics
            other_exogs: pd.DataFrame, the table encoding other traits,
                         the columns is traits, rows are individuals, values are numerics
            endog: pd.DataFrame or pd.Series, the target/endogenous variable
            engine: statsmodels models, most of the time it's
                    statsmodels.discrete.discrete_model.Logit or
                    statsmodels.regression.linear_model.OLS,
                    but any model's .fit results provide .params .bse .pvalues will work
            one_to_one_exogs: pd.DataFrame, the dataframe providing
                              specific exog variable based on encoded_snp
            one_to_one_strategy: funcs, when given a snp name, this function should output
                                 the corresponding one_to_one_exogs column name
            one_to_one_alias: list[str], the name used to take the one-to-one column
            snps_preprocessing_strategy: Optional[funcs] the preprocessing applied to snps
                   before the regression analysis
                   this function should take a pd.DataFrame as input and output
            snp_alias: str, the name used for snp column
            na_strategy: Optional[str], stands for the strategy dealing with NAs,
                 Available options are ‘none’, ‘drop’, and ‘raise’.
                 If ‘none’, no nan checking is done.
                 If ‘drop’, any observations with nans are dropped.
                 If ‘raise’, an error is raised. Default is ‘none’.
            group_na_strategy: Optional[str], how na_strategy apply when have
                              extra_iterate_on, can be "snp-wise" or "group-wise",
                              it's working for status when na_strategy == "drop"
            extra_label: Optional[tuple[tuple[str], tuple[str]]]

        Returns:
            tuple containing

            - pd.DataFrame: the effect size, se and p-values for each <endog>~column+<other> on <engine>.
            - dict: additional information that is included.
        """
        # the list saving output values
        params, bse, p_values, processed_snp, n_sample = [], [], [], [], []
        model_metric_1, model_metric_2, model_metric_3, model_metric_4 = [], [], [], []
        errors = {}
        rel_freqs, abs_freqs, totals = [], [], []
        # if there's a given na_strategy as group-wise drop NA
        # for the given group, drop every row with NA
        if (na_strategy == "drop") & (group_na_strategy == "group-wise"):
            # drop NA's
            snp_table.dropna(axis=0, how="any", inplace=True)
            other_exogs.dropna(axis=0, how="any", inplace=True)
            endog.dropna(axis=0, how="any", inplace=True)
            # when have one-to-one, drop the one-to-one as well
            if one_to_one_strategy is not None:
                one_to_one_exogs.dropna(axis=0, how="any", inplace=True)
        # Otherwise the drop in each regression(i.e, based on each SNPs' data availability)
        # for columns run the regression
        for snp in tqdm(snp_table.columns, position=0, leave=False):
            # this try is to capture regression running errors
            try:
                # prepare the exog table
                if one_to_one_strategy is None:
                    exog = pd.concat([snp_table[snp], other_exogs], join="inner", axis=1)
                else:
                    # if one-to-one, find the one-to-one column
                    one_to_one_columns = one_to_one_strategy(snp)
                    other_one_to_one_actual = one_to_one_exogs[[one_to_one_columns]]
                    exog = pd.concat(
                        [snp_table[snp], other_exogs, other_one_to_one_actual],
                        join="inner",
                        axis=1)
                # this drop is unnecessary from regression perspective,
                # but statsmodels.models cannot return nan mask now
                # so, for calculating frequencies-related, use the following
                if na_strategy == "drop":
                    exog = exog.dropna(axis=0)
                    endog_dropped = endog.dropna(axis=0)
                else:
                    endog_dropped = endog
                self.debug_log("exogs", exog, "endogs", endog_dropped)
                # align the exog and endog table
                exog, endog_dropped = exog.align(endog_dropped, join="inner", axis=0)
                # apply the preprocessing function to snps(filter_c)
                if snps_preprocessing_strategy is not None:
                    # find the snp column from current exog table
                    # apply the snps_preprocessing_strategy given
                    self.debug_log("for SNP preprocessing", exog[[snp]])
                    self.debug_log("distribution", exog[[snp]].value_counts())
                    snp_table_exog = snps_preprocessing_strategy(exog[[snp]])
                    self.debug_log("after SNP preprocessing", snp_table_exog)
                    # if no column is left
                    if len(snp_table_exog.columns) == 0:
                        raise ValueError("filtered out by snps_preprocessing_strategy")
                rel_freq, abs_freq, total = self._calculate_freq(exog, snp)
                # initialize <engine> for regression
                # if na_strategy == "drop" and group_na_strategy is "snp-wise"
                # the sample will be dropped here
                # if it's group-wise, the input should have no NAs
                # ,and it's fine to run missing = drop
                self.debug_log("exogs", exog, "endogs", endog_dropped, "corr", exog.corr())
                # This astype(float) will consume a lot of memory, but it's a must
                # otherwise statsmodels.Logit and OLS will throw an error
                # see https://stackoverflow.com/q/33833832
                regression = engine(
                    exog=exog.astype(float),
                    endog=endog_dropped.astype(float),
                    missing=na_strategy)
                # exog.to_excel(f"label_{extra_label[1]}_exog.xlsx")
                # endog.to_excel(f"label_{extra_label[1]}_endog.xlsx")
                self.debug_log("regression", regression)
                # run the regression
                result = regression.fit(disp=0)
                self.debug_log("result", result)
                self.debug_log("result summary", result.summary())
                # record the result
                params.append(result.params)
                bse.append(result.bse)
                p_values.append(result.pvalues)
                processed_snp.append(snp)
                n_sample.append(result.nobs)
                # add metrics evaluating model performances
                if isinstance(result, BinaryResultsWrapper):
                    count = pd.concat([exog[snp], endog_dropped], axis=1).value_counts()
                    model_metric_1.append(count[(1.0, 1.0)] if count.index.isin([(1.0, 1.0)]).any() else 0)
                    model_metric_2.append(count[(1.0, 0.0)] if count.index.isin([(1.0, 0.0)]).any() else 0)
                    model_metric_3.append(count[(0.0, 1.0)] if count.index.isin([(0.0, 1.0)]).any() else 0)
                    model_metric_4.append(count[(0.0, 0.0)] if count.index.isin([(0.0, 0.0)]).any() else 0)
                elif isinstance(result, RegressionResultsWrapper):
                    model_metric_1.append(result.rsquared)
                # when asking for meta info, provide them
                rel_freqs.append(rel_freq)
                abs_freqs.append(abs_freq)
                totals.append(total)
                # Memory efficiency
                gc.collect()
            except Exception as e:
                # record errors in errors[snp]
                self.debug_log(e)
                errors[snp] = str(e)
        # the name of all variables
        variable_name = exog.columns.tolist()
        variable_name[0] = snp_alias
        # when have one to one columns
        if one_to_one_strategy is not None:
            # replace the one-to-one column's name with name given
            variable_name[-other_one_to_one_actual.shape[1]:] = one_to_one_alias
        # tidy up the format
        df_result = pd.DataFrame(
            [pd.concat([params[i], bse[i], p_values[i]], axis=1).to_numpy().flatten()
             for i in range(len(params))],
            columns=self._get_column_name(variable_name),
            index=processed_snp
        )

        try:
            # only trying to capture if result exist...
            result
        except Exception as e:
            print(e)
            raise RuntimeError(
                "no valid result is recorded, usually because of all the SNPs are filtered out by "
                "snps_preprocessing_strategy")

        if isinstance(result, statsmodels.discrete.discrete_model.BinaryResultsWrapper):
            df_result["with_SNP_with_trait"] = model_metric_1
            df_result["with_SNP_without_trait"] = model_metric_2
            df_result["without_SNP_with_trait"] = model_metric_3
            df_result["without_SNP_without_trait"] = model_metric_4
            df_result["chi_sq_p_value"] = df_result.apply(
                lambda x: scipy.stats.chi2_contingency(
                    [[x["with_SNP_with_trait"], x["with_SNP_without_trait"]],
                     [x["without_SNP_with_trait"], x["without_SNP_without_trait"]]],
                    correction=False
                )[1], axis=1)
        elif isinstance(result, statsmodels.regression.linear_model.RegressionResultsWrapper):
            df_result["r_squared"] = model_metric_1
        else:
            # only left here for completeness in case need other sanity check for other result,
            # write here
            pass
        df_result["rel_freqs"] = rel_freqs
        df_result["abs_freqs"] = abs_freqs
        df_result["n_sample"] = n_sample
        if len(n_sample) == 0:
            raise ValueError("no n_sample is recorded")
        # generate an info dictionary
        analysis_info = self._generate_info(
            endog=endog,
            engine=engine.__name__,
            variable_name=variable_name,
            extra_label=extra_label,
            processed_snp=processed_snp,
            n_sample=n_sample[-1],
            errors=errors)
        return df_result, analysis_info

    def association_test(
            self,
            encoded_snp: pd.DataFrame,
            other_exogs: pd.DataFrame,
            target_dataset: Union[pd.DataFrame, pd.Series],
            target_strategy: Dict[str, Dict[str, Union[sm.Logit, sm.OLS, Callable]]],
            output_path: str = "",
            one_to_one_exogs: Optional[pd.DataFrame] = None,
            one_to_one_strategy: Optional[Callable] = None,
            one_to_one_alias: str = None,
            na_strategy: str = "drop",
            extra_iterate_on: List[str] = None,
            group_na_strategy: str = "snp_wise",
            snps_preprocessing_strategy: Optional[Callable] = None,
            snp_alias: str = "variant",
            verbose: int = 0
    ):
        """API running regression test

        Args:

            encoded_snp: pd.DataFrame, the dataframe to be looped on columns
            other_exogs: pd.DataFrame, the dataframe taking all the other variables
            target_dataset: pd.DataFrame, the dataframe taking all the target variables
            target_strategy: dict[str, dict[str, funcs or models]], the dictionary
                provide pre-processing to specific column. Only column mentioned in keys
                will be included in the running. The inner dictionary should have two keys:

                 - "engine": statsmodels.api models,
                     designed for statsmodels.discrete.discrete_model.Logit or
                     statsmodels.regression.linear_model.OLS,
                     but any model's .fit results provide .params .bse .pvalues will work

                 - "preprocessing": funcs
                     the function should take a pd.DataFrame/pd.Series as input
                     and a pd.DataFrame/pd.Series as output,
                     This strategy is None, the column will be used as
                     exogenous variables as-is

            output_path: str, the output root path, the saving name will be printed for reference
            one_to_one_exogs: pd.DataFrame, the dataframe providing
                            specific exog variable based on encoded_snp
            one_to_one_strategy: funcs, when given a snp name, this function should output
                                 the corresponding one_to_one_exogs column name
            one_to_one_alias: str, the name for one_to_one settings
            snps_preprocessing_strategy: Optional[funcs] the preprocessing applied to snps
                           before the regression analysis
                           this function should take a pd.DataFrame as input and output
            na_strategy: Optional[str], stands for the strategy dealing with NAs,
                 Available options are ‘none’, ‘drop’, and ‘raise’.
                 If ‘none’, no nan checking is done.
                 If ‘drop’, any observations with nans are dropped.
                 If ‘raise’, an error is raised. Default is ‘none’.
            extra_iterate_on: Optional[list], the list of extra iterate variable in other_exog
                              please be noticed that it may cause severe explosions in running time
            group_na_strategy: Optional[str], how na_strategy apply when have
                               extra_iterate_on, can be "snp-wise" or "group-wise",
                               it's working for status when na_strategy == "drop"
            snp_alias: Optional[str], the name used for snp column
            verbose: Optional[int] in {0,1,2} default 0
                     if verbose = 1, the regression will give the saving path
                     if verbose = 2, the regression will output all the related
                     values during the computation, only for debugging purpose,
                     use with care for it will create massive I/O
        Returns:
            dict, a dict recording output brief information
        """
        os.makedirs(output_path, exist_ok=True)
        output_info = {}
        # iterate over all the variables mentioned
        for target, strategy in tqdm(target_strategy.items(), position=0, leave=True):
            output_info[target] = {}
            # refresh the memory
            gc.collect()
            # decode the strategy given
            preprocessing = strategy["preprocessing"]
            # if there's no strategy given
            if preprocessing is None:
                # use the original column
                endogenous = target_dataset[[target]]
            # otherwise apply strategy function given
            else:
                endogenous = preprocessing(target_dataset[[target]])
            # when not given extra iterates
            if extra_iterate_on is None:
                # run the association analysis
                result, info = self._association_snp(
                    snp_table=encoded_snp,
                    other_exogs=other_exogs,
                    endog=endogenous,
                    engine=strategy["engine"],
                    one_to_one_exogs=one_to_one_exogs,
                    one_to_one_strategy=one_to_one_strategy,
                    one_to_one_alias=one_to_one_alias,
                    snps_preprocessing_strategy=snps_preprocessing_strategy,
                    snp_alias=snp_alias,
                    na_strategy=na_strategy,
                    group_na_strategy=group_na_strategy,
                    extra_label=None
                )
                self._save_association_result(result, info, output_path, verbose)
                output_info[target]["output"] = info
            # when given extra iterates
            else:
                # check all the extra_iterates
                other_exogs_group = other_exogs.groupby(extra_iterate_on)
                for sub_exogs in tqdm(other_exogs_group, position=1, leave=False):
                    result, info = self._association_snp(
                        snp_table=encoded_snp,
                        other_exogs=sub_exogs[1].drop(columns=extra_iterate_on),
                        endog=endogenous,
                        engine=strategy["engine"],
                        one_to_one_exogs=one_to_one_exogs,
                        one_to_one_strategy=one_to_one_strategy,
                        one_to_one_alias=one_to_one_alias,
                        snps_preprocessing_strategy=snps_preprocessing_strategy,
                        snp_alias=snp_alias,
                        na_strategy=na_strategy,
                        group_na_strategy=group_na_strategy,
                        extra_label=(extra_iterate_on, [sub_exogs[0]])
                    )
                    self._save_association_result(result, info, output_path, verbose)
                    output_info[target][sub_exogs[0]] = info
        gc.collect()
        return output_info
