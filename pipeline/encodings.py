"""All the settings related to encodings"""
import pandas as pd
import os
import gc
import glob
import warnings

def is_notebook() -> bool:
    """https://stackoverflow.com/a/39662359"""
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

if is_notebook():
    from tqdm.notebook import tqdm
else:
    from tqdm.auto import tqdm
tqdm.pandas()


class EncodingCoassinOutput():
    """
    The pipeline running encoding outputs from coassin_pipeline

    Encoding Rule:
        If Raw file total coverage <50, make that position missing for that person
        If Raw file total coverage >=50, then look at annotated file
            If position is missing in annotated file, the variant is coded 0
            If position is present in the annotated, variant_level>0.01 and
            the total reads supporting (variant_level*total_coverage) the variant
            be >= 10,the variant is coded 1, otherwise 0

    Algorithm:
    1. Read and encode each folders individually,
       if a snp position didn't give mutation result, save them as pos, else pos-ref/alt,
    2. When we get all the encodings, concate them together and deal with all the pos conflict

    Hardware Requirement:
    The algorithm is preventing the growth of the number of cross-table accesses
    and queries, especially minimizing the time read in dataframe from the file
    system(= 2*n_samples) and its memory cost (a raw.txt and a variantsAnnotate.txt
    at any specific time) during 1st step encoding. Which is at the cost of
    HUGE memory requirement dealing with concatenation (which is more customizable
    on our cluster).
    For reproducing, generate_encoded_results() method need at least 10GB memories.
    For a smooth running, I suggest at least 20GB to 30GB memories for the kernel.

    Initialize:
        eco = EncodingCoassinOutput(choose 1: input_path = "/some/parent/path/of/bam/output" or
                                              bam_list = "/paths/to/a/file/recording/bam/path/line/by/line.txt"
                                    output_path = "output_path", # required
                                    raw_total_coverage_threshold = 50,
                                    variant_level_threshold = 0.01,
                                    supporting_threshold = 10)

    Encoding:
        1. run eco.encode_individual(saving_step = 500(or any integer)),
        this step will be very time consuming. For it is originally running on SGE,
        no parallel is provided in Python, it's recommended to split your BAM output via bam_list
        After 1 finish, run the following in one thread:
        2. eco.generate_coverage_total() for final coverage_total
        3. eco.generate_encoded_results() for final encoded_results
    """
    def __init__(self,
                 output_path,
                 input_path = None,
                 bam_list = None,
                 raw_total_coverage_threshold = 50,
                 variant_level_threshold = 0.01,
                 supporting_threshold = 10,
                 verbose = True):
        """
        The initializer inputing all the results
        Args:
            output_path: str, the path saving all encoding output
            input_path: Optional[str], the path saving all the coassin's output
            bam_list: Optional[list], the name of all the coassin output we will use
            raw_total_coverage_threshold: int, see above.
            variant_level_threshold: float, see above.
            supporting_threshold: int, see above.
            verbose: bool, default True: if False, no log will be printed
        """
        # save the output_path, create one if not provided
        self._output_path = output_path
        os.makedirs(self._output_path, exist_ok = True)
        # if the bam_list is not provided, search over the folders input_path
        if bam_list is None:
            # glob.iglob("input_path/**/*.BQSR.recaled.bam", recursive = True)
            # searches every file and subfolder under input_path
            bam_list = glob.iglob(os.path.join(
                                    input_path,
                                    "**",
                                    "*.BQSR.recaled.bam"),
                                    recursive = True)
        else:
            with open(bam_list, "r") as file:
                bam_list = [line.rstrip() for line in file]
        self._verbosity = verbose
        # Verify that each folder provided have variantsAnnotate/variantsAnnotate.txt file
        self._bam_list = [bam_output for bam_output in tqdm(bam_list)
                          if self._verify_coassin_output(bam_output)]
        if self._verbosity == True:
            print(f"{len(self._bam_list)} valid input detected")
        del bam_list
        gc.collect()

        self._raw_total_coverage_threshold = raw_total_coverage_threshold
        self._variant_level_threshold = variant_level_threshold
        self._supporting_threshold = supporting_threshold
        pd.set_option('mode.chained_assignment', None)
        # this warning is triggered by generate_encoded_results() method, and is properly dealt with
        warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


    def _verify_coassin_output(self, path):
        """
        assert <path> is a complete Coassin path output

        check the existence of
        1. <path>/raw/raw.txt
        2. <path>/variantsAnnotate/variantsAnnotate.txt

        Args:
            path: str, the path to a Coassin pipeline output
        Returns
            bool, if the output is valid, return True, otherwise False
        """
        #
        return (os.path.isfile(
          os.path.join(path,"raw","raw.txt")
        ) and os.path.isfile(
          os.path.join(path,"variantsAnnotate","variantsAnnotate.txt")
        ))

    def _encode_individual(self, bam_path):
        """The procedure reading in, encoding, then tidy up the output for an invidual"""
        # read the coverage total column from raw.txt file
        coverage_total, SampleID = self._get_coverage_total(
            self._data_cleaning_raw(
                self._get_file_raw(bam_path = bam_path)))
        # read the variant_annotated.txt file
        annotated_files = self._data_cleaning_annotated(
            self._get_file_annotated(bam_path = bam_path))
        # group the coverage_total by its position for row-wise encoding
        encoded_group = coverage_total.reset_index().groupby(["POS","COV-TOTAL"])
        encode_result = encoded_group.apply(
                            lambda x: self._encode_position(
                                pos = x.name[0],
                                coverage_total = x.name[1],
                                annotated_files = annotated_files))
        # tidy up for merge and output
        coverage_total = coverage_total.rename(columns = {"COV-TOTAL": SampleID})
        encode_result = encode_result.reset_index().drop(
                            columns = ["COV-TOTAL","POS","level_2"]
                            ).set_index("snp_pos").astype(pd.Float64Dtype())
        encode_result = encode_result.rename(columns = {"encoding": SampleID})
        # clean the memory
        gc.collect()
        return coverage_total, encode_result

    def _encode_position(self, pos, coverage_total, annotated_files):
        """actual encoding logic, apply to each specific row(allele) in variantsAnnotate file

        both pos and coverage_total is obtained from pd.DataFrameGroupBy object's apply method

        Args:
            pos: int, the position of this SNP
            coverage_total: int, the coverage_total read from the raw.txt
            annotated_files: pd.DataFrame, reading from
                             variantsAnnotate/variantsAnnotate.txt
        """

        # if the coverage_total is less than raw_total_coverage_threshold
        if coverage_total < self._raw_total_coverage_threshold:
            # encode it as pos, pd.NA
            return pd.DataFrame([[pos, pd.NA]], columns = ["snp_pos", "encoding"])
        else:
            #If position is present in the annotated, then variant_level>0.01
            # and the total reads supporting (variant_level*total_coverage) the variant should be>=10
            if pos in annotated_files["Pos"].values:
                related_annotation = annotated_files[annotated_files["Pos"].eq(pos)]
                # variant_level>0.01 and the total reads supporting
                # (variant_level*total_coverage) the variant should be>=10
                related_annotation["encoding"] = \
                    (related_annotation["Variant-Level"]>self._variant_level_threshold) & \
                    (related_annotation["Variant-Level"]*related_annotation["Coverage-Total"]\
                     >self._supporting_threshold)
                # encode the True or False to 1 or 0
                return related_annotation[["snp_pos", "encoding"]].replace({True: 1, False: 0})
            # if the annotated files doesn't have this row, encode them as 0
            else:
                return pd.DataFrame([[pos, 0]], columns = ["snp_pos", "encoding"])

    def _tidy_encoded_results(self, df):
        """read in a 1st step encoded_result from path, then tidy the format up

        inplace=True is applied to minimize the memory cost
        Args:
            path: str, the path of the file
        Returns:
            df: pd.DataFrame, the encoded_result
        """
        df.fillna(-1, inplace = True)
        df.reset_index(inplace = True)
        df.drop_duplicates(inplace = True)
        df.set_index("snp_pos", inplace = True)
        return df

    def encode_individual(self, saving_step = 500):
        """API for the first step individual-wise encoding procedure,

        It will take very long time, so for loop and huge amount of intermediate
        results will be saved in self._output_path, to customize the saving step,
        a saving_step (default 500) is provided

        Args:
            save_step: int, after every <save_step> sample, the program will save the result
        Saves:
            self._output_path/coverage_totals/<sample_name>.csv
            self._output_path/encoded_results/<sample_name>.csv
        """
        # create the saving path
        coverage_totals_path = os.path.join(self._output_path, "coverage_totals")
        os.makedirs(coverage_totals_path, exist_ok=True)
        encode_results_path = os.path.join(self._output_path, "encoded_results")
        os.makedirs(encode_results_path, exist_ok=True)
        sample_size = len(self._bam_list)
        # loop through the bam_list
        for saving_loop in tqdm(range(0, sample_size, saving_step), position=0, leave=True):
            # encode the slice of bam list, save the result into two individual lists
            coverage_total_list, encode_result_list = \
                zip(*(self._encode_individual(bam_path = self._bam_list[saving_loop*saving_step + i])
                    for i in range(min(saving_step, sample_size - saving_loop*saving_step))))
            # save coverage total and encoding result
            coverage_total_combined = pd.concat(coverage_total_list, axis = 1, join = "outer")
            coverage_total_saving_path = os.path.join(coverage_totals_path,
                                                      f"{os.path.basename(self._bam_list[saving_loop])}.csv")
            coverage_total_combined.to_csv(coverage_total_saving_path)
            if self._verbosity == True:
                print(f"Step: {int(saving_loop/saving_step)}: coverage total saved into {coverage_total_saving_path}")
            del coverage_total_list, coverage_total_combined
            gc.collect()
            encode_result_combined = pd.concat([self._tidy_encoded_results(encode_result)
                                                for encode_result in encode_result_list],
                                                axis = 1, join = "outer")
            encoded_result_saving_path = os.path.join(encode_results_path,
                                                      f"{os.path.basename(self._bam_list[saving_loop])}.csv")
            encode_result_combined.to_csv(encoded_result_saving_path)
            if self._verbosity == True:
                print(f"Step: {int(saving_loop/saving_step)}: encoded result saved into {encoded_result_saving_path}")
            # memory cleaning
            del encode_result_list, encode_result_combined
            gc.collect()

    def generate_coverage_total(self, save = False):
        """API generate final coverage_total table
        Args:
            save: Bool, default False, if True, it will save the final result at
                  self._output_path/coverage_total_final.csv
        Returns:
            coverage_total: pandas.DataFrame, the final coverage_total table
        Saves:
            if save == True, self._output_path/coverage_total_final.csv, same as above
        """
        ct_path_iter = glob.iglob(os.path.join(self._output_path, "coverage_totals", "**","*.csv"),
                                  recursive = True)
        ct_iter = (pd.read_csv(path, index_col = 0) for path in ct_path_iter)
        coverage_total = pd.concat(ct_iter, axis = 1, join = "outer")
        if save == True:
            coverage_total.to_csv(os.path.join(self._output_path, "coverage_total_final.csv"))
        # Memory Efficiency
        del ct_path_iter, ct_iter
        gc.collect()
        return coverage_total

    def generate_encoded_results(self, save = False, tidy_when_load = False):
        """The procedure generate final encoded_results table

        Args:
            save: Bool, default False, if True, it will save the final result at
                  self._output_path/encoded_result_final.csv
            tidy_when_load: Bool, default False, if True, when load each seperate encoded_result
                            will apply self._tidy_encoded_results, just for running different version code
                            please use False in your application.
        Returns:
            combined_er: pandas.DataFrame, the final encoded_result table
        Saves:
            if save == True, self._output_path/encoded_result_final.csv, same as above
        """
        # a list of all the encoded result
        er_path_iter = glob.iglob(os.path.join(self._output_path,"encoded_results", "**","*.csv"),
                                  recursive = True)
        # for each encoded results tidy up the format
        if tidy_when_load == True:
            er_iter = (self._tidy_encoded_results(pd.read_csv(path, index_col = 0))
                       for path in tqdm(er_path_iter, position=0, leave=True))
        else:
            er_iter = (pd.read_csv(path, index_col = 0)
                       for path in tqdm(er_path_iter, position=0, leave=True))
        # concate the encoded result
        er_complete = pd.concat(er_iter, join = "outer", axis = 1)
        # deleted the huge list with csvs
        del er_path_iter, er_iter
        gc.collect()
        # generate a "position" column
        er_complete["pos"] = er_complete.index.map(lambda x: str(x).split("-")[0])
        # please keep this copy is
        # 1. to combine the dataframe as a continuous chunk in RAM
        # 2. necessary for modifying the Dtypes,
        er_complete = er_complete.copy()
        gc.collect()
        er_complete.index = er_complete.index.astype(pd.StringDtype())
        er_complete["pos"] = er_complete["pos"].astype(pd.StringDtype())
        # split the whole table by the pos
        er_complete_group = er_complete.groupby("pos")
        # apply the combine wrapper to each position
        combined_er = er_complete_group.progress_apply(self._combine_position_wrapper)
        # clean the er_complete_group
        del er_complete_group
        gc.collect()
        # some formatting
        combined_er = combined_er.reset_index()
        # sort the output by position
        combined_er["pos"] = combined_er["pos"].apply(lambda x: int(x))
        combined_er = combined_er.sort_values("pos")
        combined_er = combined_er.drop(columns = ["index", "pos"], errors = "ignore")
        combined_er.set_index("snp_pos", inplace = True)
        if save == True:
            combined_er.to_csv(os.path.join(self._output_path, "encoded_result_final.csv"))
        return combined_er

    def _combine_position_wrapper(self, x):
        """wrapper function apply self._combine_encoded_on_position() to a DataFrame"""
        base = x.name
        x = x.apply(self._combine_encoded_on_position, axis = 0, base = base)
        x.drop(index = [f"{base}", base], columns = "pos",inplace= True,errors = "ignore")
        gc.collect()
        return x

    def _combine_encoded_on_position(self, col, base):
        """the actual logic combine the result of a column

        Run the following logic:

        if "position" appears
        case 1:
            if "position" is "missing", total coverage < 50, the final encoding should be NA
        case 2:
            if both "position" and all "position-ref/alt" is NA:
            raw.txt doesn't provide this position, the final encoding should be all NA
        case 3:
            if "position" is encoded, and "position-ref/alt" is NA
            raw.txt provided this value, total coverage > 50,
            but this "position-ref/alt" is not provided in annotatedVariant.txt
            The final encoding should be the value in "position" = 0
            only "position-ref/alt" can be encoded as 1
        case 4:
            if "position" is NA, and "position-ref/alt" is encoded
            raw.txt provided this value
            and this "position-ref/alt" is provided in annotatedVariant.txt
            The final encoding should be the value in "position-ref/alt",
            if "position-ref/alt" is NA but with other Non-NA "position-ref/alt",
            it's 0 for any NA for the total coverage is proven enough
        case 5:
            Both "position" and "position-ref/alt" are encoded:
            Its technically not possible
        if "position" is not there, only "position-ref/alt" appears
        case 6:
            If any of these "position-ref/alt" is encoded:
            Encoded position-ref/alt keeps the value in "position-ref/alt"
            All the rest position-ref/alt will be 0, they are not provided in annotatedVariant.txt,
            but the coverage total is proven > 50
        case 7:
            If all "position-ref/alt" is NA: The final encoding should be 0,
            the coverage total is proven > 50
        """
        # convert base to str, just format problem
        base = f"{base}"
        # if a "<position>" is in the case, "position" will only be 0/missing/NA
        if base in col.index:
            # find the row whose index is the "base", i.e. "<position>",
            # rather than <position>-<ref>/<alt>
            base = col.loc[base]
            # if coverage_total < 50 or all NA
            if (base == "missing") or (col.isna().all()) or (base == -1):
                # Case 1: return as-is, base will be dropped, all the rest is natually NA
                # Case 2: all the element is natually NA
                return col
            # if coverage > 50 and not all NA, base = 0 or NA
            else:
                # Case 3: if position is encoded 0, raw.txt provided a total coverage > 50
                # but not in annotatedVariant.txt all the rest will be NA, so fill them as 0
                # Case 4: otherwise position has to be NA, and the col has some non-NA position-ref/alt
                return col.fillna(0)
        # if there's no "position", total coverage > 50 for sure
        else:
            #case 6 and 7
            return col.fillna(0)


#----------------------------reading files--------------------------------------
    def _extract_ID(self, SampleID):
        """
        Helper function clean Sample ID to pure digits

        Args:
            SampleID: String, in "washei*****.BQSR.recaled.bam" format, where * stand for numbers
        Return:
            String, the WES ID (the part before ".BQSR..." )
        """
        return SampleID.split(".")[0]

    def _get_file_annotated(self, bam_path):
        '''
        given a bam name, read the variantsAnnotate.txt inside the output of Coassin pipeline
        Args:
            input_path, the input path with all the bam output inside
            bam_path: str, the name of the original bam file
        Return:
            pandas.DataFrame instance, the variantsAnnotate.txt file read
        '''
        return pd.read_csv(os.path.join(bam_path,"variantsAnnotate","variantsAnnotate.txt"), delimiter = "\t")

    def _data_cleaning_annotated(self,df):
        """
        existing variantsAnnotate.txt data cleaning procedure
        Args:
            df: pandas.DafaFrame instance, the dataframe
        Return:
            df: pandas.DafaFrame instance, the dataframe with cleaned ID
        """
        df["SampleID"] = df["SampleID"].apply(self._extract_ID)
        # prepare a snp_pos format label("pos-ref/alt") for annotated_files
        df["snp_pos"] = df.apply(lambda x: f"{x['Pos']}-{x['Ref']}/{x['Variant']}", axis = 1)
        return df

    def _get_file_raw(self, bam_path):
        '''
        given a bam name, read the <bam_path>/raw/raw.txt inside the output of Coassin pipeline
        Args:
            input_path, the input path with all the bam output inside
            bam_path: str, the name of the original bam file
        Return:
            pandas.DataFrame instance, the variantsAnnotate.txt file read
        '''
        return pd.read_csv(os.path.join(bam_path,"raw","raw.txt"), delimiter = "\t")

    def _data_cleaning_raw(self, df):
        """
        existing raw.txt data cleaning procedure
        Args:
            df: pandas.DafaFrame instance, the dataframe
        Return:
            df: pandas.DafaFrame instance, the dataframe with cleaned ID
        """
        df["SAMPLE"] = df["SAMPLE"].apply(self._extract_ID)
        df = df.rename(columns = {"SAMPLE": "SampleID"})
        return df

    def _get_coverage_total(self, raw):
        # get the Sample ID
        SampleID = raw.SampleID.loc[0]
        # Pick POS and coverage total, set index as POS, rename the coverage total as SampleID
        raw = raw[["POS", "COV-TOTAL"]]
        raw = raw.set_index("POS")
        raw = raw.sort_index()
        return raw, SampleID
