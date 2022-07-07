"""All the settings related to encodings"""
import numpy as np
import pandas as pd
import os
import gc
import statsmodels.api as sm
import datetime
import sys
import glob
import shutil
import itertools
from tqdm import tqdm

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
        eco = EncodingCoassinOutput(input_path = "some_path",
                                    output_path = "other_path",
                                    raw_total_coverage_threshold = 50,
                                    variant_level_threshold = 0.01,
                                    supporting_threshold = 10)
    Encoding:
        After eco.encode_individual() (ETA ~6hours)
        Run eco.generate_coverage_total() for final coverage_total
        eco.generate_encoded_results() for final encoded_results
    """
    def __init__(self,
                 input_path,
                 output_path,
                 bam_list = None,
                 raw_total_coverage_threshold = 50,
                 variant_level_threshold = 0.01,
                 supporting_threshold = 10):
        """
        The initializer inputing all the results
        Args:
            input_path: str, the path saving all the coassin's output
            output_path: str, the path saving all encoding output
            bam_list: Optional[list], the name of all the coassin output we will use
            raw_total_coverage_threshold: int, see above.
            variant_level_threshold: float, see above.
            supporting_threshold: int, see above.
        """
        self._input_path = input_path
        self._output_path = output_path
        os.makedirs(self._output_path, exist_ok = True)
        # if the bam_list is provided
        if not bam_list is None:
            # use this bam_list
            self._bam_list = bam_list
        # otherwise search over the folders input_path and its subfolder called ./hispanic
        else:
            # get its sub-folder
            input_path_hispanic = os.path.join(self._input_path, "hispanic")
            # find all the folder with proper name in input_path
            bam_list_0 = [dir for dir in next(os.walk(self._input_path))[1]
                          if '.BQSR.recaled.bam' in dir]
            # find all the folder with proper name in input_path/hispanic
            bam_list_1 = [os.path.join("hispanic",dir)
                          for dir in next(os.walk(input_path_hispanic))[1]
                          if '.BQSR.recaled.bam' in dir]
            # combine the list
            bam_list = bam_list_0+bam_list_1
            self._bam_list = [i for i in bam_list
                              if os.path.isfile(
                                os.path.join(input_path, i,
                                            "variantsAnnotate",
                                            "variantsAnnotate.txt")
                              )]
            del bam_list, bam_list_0, bam_list_1
            gc.collect()

        self._raw_total_coverage_threshold = raw_total_coverage_threshold
        self._variant_level_threshold = variant_level_threshold
        self._supporting_threshold = supporting_threshold
        pd.set_option('mode.chained_assignment', None)

    def _encode_individual(self, input_path, bam_name):
        """The procedure reading in, encoding, then tidy up the output for an invidual"""
        # read the coverage total column from raw.txt file
        coverage_total, SampleID = self._get_coverage_total(
            self._data_cleaning_raw(
                self._get_file_raw(
                    input_path = input_path,
                    bam_name = bam_name)))
        # read the variant_annotated.txt file
        annotated_files = self._data_cleaning_annotated(
            self._get_file_annotated(
                input_path = input_path,
                bam_name = bam_name))
        # prepare a snp_pos format label for annotated_files
        annotated_files["snp_pos"] = annotated_files.apply(lambda x: f"{x['Pos']}-{x['Ref']}/{x['Variant']}", axis = 1)
        # group the coverage_total by its position for row-wise encoding
        encoded_group = coverage_total.reset_index().groupby("POS")
        # apply the encoding procedure
        encode_result = encoded_group.apply(self._encode_position, annotated_files)

        # tidy up for merge and output
        coverage_total = coverage_total.rename(columns = {"COV-TOTAL": SampleID})
        encode_result = encode_result.reset_index().drop(columns = ["POS","level_1"]).set_index("snp_pos")
        encode_result = encode_result.rename(columns = {"encoding": SampleID})
        # clean the memory
        gc.collect()
        return coverage_total, encode_result

    def _encode_position(self, row, annotated_files):
        """actual encoding logic, apply to each specific row(allele) in variantsAnnotate file

        Args:
            row: pd.DataFrame, a dataframe on a specific pos,
                 filtered by pd.groupby("POS")
            annotated_files: pd.DataFrame, reading from
                             variantsAnnotate/variantsAnnotate.txt
        """
        # get the position and coverage_total
        pos = row["POS"].iloc[0]
        coverage_total = row["COV-TOTAL"].iloc[0]
        # if the coverage_total is less than raw_total_coverage_threshold
        if coverage_total < self._raw_total_coverage_threshold:
            # encode it as pos, pd.NA
            return pd.DataFrame([[pos, pd.NA]], columns = ["snp_pos", "encoding"])
        else:
            #If position is present in the annotated, then variant_level>0.01
            # and the total reads supporting (variant_level*total_coverage) the variant should be>=10
            if sum(annotated_files["Pos"].eq(pos)):
                related_anno = annotated_files[annotated_files["Pos"].eq(pos)]
                # variant_level>0.01 and the total reads supporting
                # (variant_level*total_coverage) the variant should be>=10
                related_anno["encoding"] = \
                    (related_anno["Variant-Level"]>self._variant_level_threshold) & \
                    (related_anno["Variant-Level"]*related_anno["Coverage-Total"]>self._supporting_threshold)
                # encode the True or False to 1 or 0
                return related_anno[["snp_pos", "encoding"]].replace([True, False], [1,0])
            # if the annotated files doesn't have this row, encode them as 0
            else:
                return pd.DataFrame([[pos, 0]], columns = ["snp_pos", "encoding"])

    def encode_individual(self):
        """API for the first step individual-wise encoding procedure,

        It will take very long time, so for loop and huge amount of intermediate
        results will be saved in self._output_path

        Saves:
            self._output_path/coverage_totals/<sample_name>.csv
            self._output_path/encoded_results/<sample_name>.csv
        """
        # create the saving path
        converage_totals_path = os.path.join(self._output_path, "coverage_totals")
        os.makedirs(converage_totals_path, exist_ok=true)
        encode_results_path = os.path.join(self._output_path, "encoded_results")
        os.makedirs(encode_results_path, exist_ok=true)
        # loop through the bam_list
        for i in tqdm(range(len(self._bam_list)),position=0, leave=True):
            # encode one person
            coverage_total, encode_result = self._encode_individual(input_path = input_path,bam_name = bam_list[i])
            # save coverage total and encoding result
            coverage_total.to_csv(os.path.join(converage_totals_path, f"{bam_list[i]}.csv"))
            encode_result.to_csv(os.path.join(encode_results_path, f"{bam_list[i]}.csv"))
            # memory cleaning
            del coverage_total, encode_result
            gc.collect()

    def generate_coverage_total(self):
        """API generate final coverage_total table

        Returns:
            coverage_total: pandas.DataFrame, the final coverage_total table
        Saves:
            self._output_path/coverage_total_final.csv, same as above
        """
        ct_path_list = glob.glob(os.path.join(self._output_path,"coverage_totals","*.csv")) + \
               glob.glob(os.path.join(self._output_path,"coverage_totals","**","*.csv"))
        ct_list = [pd.read_csv(path, index_col = 0)
                   for path in tqdm(ct_path_list,position=0, leave=True)]
        coverage_total = pd.concat(ct_list, axis = 1, join = "outer")
        coverage_total.to_csv(os.path.join(self._output_path, "coverage_total_final.csv"))
        # Memory Efficiency
        del ct_name_list, ct_list, coverage_total
        gc.collect()
        return coverage_total

    def generate_encoded_results(self):
        """The procedure generate final encoded_results table

        Returns:
            combined_er: pandas.DataFrame, the final encoded_result table
        Saves:
            self._output_path/encoded_result_final.csv, same as above
        """
        # a list of all the encoded result
        er_path_list = glob.glob(os.path.join(self._output_path,"encoded_results","*.csv")) + \
               glob.glob(os.path.join(self._output_path,"encoded_results","**","*.csv"))
        # for each encoded results tidy up the format
        er_list_cleaned = [self._tidy_encoded_results(path = path)
                           for path in tqdm(er_path_list, position=0, leave=True)]
        gc.collect()
        # concate the encoded result
        er_complete = pd.concat(er_list_cleaned, join = "outer", axis = 1)
        # deleted the huge list with csvs
        del er_list_cleaned
        gc.collect()
        # generate a "pos" column
        er_complete.reset_index(inplace=True)
        er_complete["pos"] = er_complete["snp_pos"].apply(lambda x: x.split("-")[0])
        er_complete.set_index("snp_pos", inplace = True)
        # this copy is necessary for modifying the Dtypes, please keep
        er_complete = er_complete.copy()
        er_complete.index = er_complete.index.astype(pd.StringDtype())
        er_complete.pos = er_complete.pos.astype(pd.StringDtype())
        # split the whole table by the pos
        er_complete_group = er_complete.groupby("pos")
        # apply the combine wrapper to each position
        combined_er = er_complete_group.apply(self._combine_position_wrapper)
        # clean the er_complete_group
        del er_complete_group
        gc.collect()
        # some formatting
        combined_er = combined_er.reset_index()
        # sort the output by position
        combined_er["pos"] = combined_er["pos"].apply(lambda x: int(x))
        combined_er = combined_er.sort_values("pos")
        combined_er = combined_er.drop(columns = ["index", "pos"], errors = "ignore")
        combined_er.set_index("snp_pos")
        combined_er.to_csv(os.path.join(self._output_path, "encoded_result_final.csv"))
        return combined_er

    def _tidy_encoded_results(self, path):
        """read in a 1st step encoded_result from path, then tidy the format up

        inplace=True is applied to minimize the memory cost
        Args:
            path: str, the path of the file
        Returns:
            df: pd.DataFrame, the encoded_result
        """
        df = pd.read_csv(path, index_col = 0)
        df.fillna("missing", inplace = True)
        df.reset_index(inplace = True)
        df.drop_duplicates(inplace = True)
        df.set_index("snp_pos", inplace = True)
        return df

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
            if (base == "missing") or (col.isna().all()):
                # Case 1: return as-is, base will be dropped, all the rest is natually NA
                # Case 2: all the element is natually NA
                return col
            # if coverage > 50 and not all NA, base = 0 or NA
            else:
                # Case 3: if position is encoded 0, raw.txt provided a total coverage > 50 but not in annotatedVariant.txt
                # all the rest will be NA, so fill them as 0
                # Case 4: otherwise position has to be NA, and the col has some non-NA position-ref/alt
                return col.fillna(0)
        # if there's no "position", total coverage > 50 for sure
        else:
            #case 6 and 7
            return col.fillna(0)

    def _combine_position_wrapper(self, x):
        """wrapper function apply self._combine_encoded_on_position() to a DataFrame"""
        base = x.name
        x = x.apply(self._combine_encoded_on_position, axis = 0, base = base)
        x.drop(index = [f"{base}", base], columns = "pos",inplace= True,errors = "ignore")
        gc.collect()
        return x
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

    def _get_file_annotated(self.input_path, bam_name):
        '''
        given a bam name, read the variantsAnnotate.txt inside the output of Coassin pipeline
        Args:
            input_path, the input path with all the bam output inside
            bam_name: str, the name of the original bam file
        Return:
            pandas.DataFrame instance, the variantsAnnotate.txt file read
        '''
        return pd.read_csv(f'{input_path}/{bam_name}/variantsAnnotate/variantsAnnotate.txt', delimiter = "\t")

    def _get_file_raw(self, input_path, bam_name):
        '''
        given a bam name, read the variantsAnnotate.txt inside the output of Coassin pipeline
        Args:
            input_path, the input path with all the bam output inside
            bam_name: str, the name of the original bam file
        Return:
            pandas.DataFrame instance, the variantsAnnotate.txt file read
        '''
        return pd.read_csv(f'{input_path}/{bam_name}/raw/raw.txt', delimiter = "\t")

    def _data_cleaning_annotated(self,df):
        """
        existing annotatedVariant data cleaning procedure
        Args:
            df: pandas.DafaFrame instance, the dataframe
        Return:
            df: pandas.DafaFrame instance, the dataframe with cleaned ID
        """
        df.SampleID = df.SampleID.apply(self._extract_ID)
        return df

    def _data_cleaning_raw(self,df):
        """
        existing raw data cleaning procedure
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
