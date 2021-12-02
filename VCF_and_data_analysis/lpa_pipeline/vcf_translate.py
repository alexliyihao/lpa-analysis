import pandas as pd
import os
import datetime
import argparse
import sys
import re

def write(path, content):
    """
    a wrapper function for with-open-write, just for good code looking...
    Args:
        path: String, the path
        content: String, the content to be written
    """
    with open(path, "w+") as f:
        read_data = f.read()
        f.write(content)

def append(path, content):
    """
    a wrapper function for with-open-append, just for good code looking...
    Args:
        path: String, the path
        content: String, the content to be written
    """
    with open(path, "a+") as f:
        read_data = f.read()
        f.write(content)

def extract_ID(SampleID):
    """
    Helper function clean Sample ID to pure digits

    Args:
        SampleID: String, in "washei*****.BQSR.recaled.bam" format, where * stand for numbers
    Return:
        String, the WES ID (the part before ".BQSR..." )
    """
    return SampleID.split(".")[0]

def get_file(input_path, bam_name):
    '''
    given a bam name, read the variantsAnnotate.txt inside the output of Coassin pipeline
    Args:
        input_path, the input path with all the bam output inside
        bam_name: str, the name of the original bam file
    Return:
        pandas.DataFrame instance, the variantsAnnotate.txt file read
    '''
    return pd.read_csv(f'{input_path}/{bam_name}/variantsAnnotate/variantsAnnotate.txt', delimiter = "\t")

def data_cleaning(df):
    """
    existing data cleaning procedure
    Args:
        df: pandas.DafaFrame instance, the dataframe
    Return:
        df: pandas.DafaFrame instance, the dataframe with cleaned ID
    """
    df.SampleID = df.SampleID.apply(extract_ID)
    return df

def data_load_wrapper(input_path, bam_list):
    """
    The wrapper for a loading files
    Args:
        path_list: list[str], the path hold all Coassin pipeline output
    Return:
        df: pandas.DafaFrame instance, all the data should be loaded
    """
    return pd.concat([
        data_cleaning(
            get_file(input_path = input_path,
                     bam_name = bam_name)
        )
        for bam_name in bam_list
    ],
        axis = "index",
        ignore_index = True)

def generate_meta(file_format = "VCFv4.2",
                  file_date = 'today',
                  source = "lpa-analysis",
                  reference = "https://raw.githubusercontent.com/seppinho/mutation-server/76e865ece25cf792d1534b0288b2c28bc1b3d013/test-data/dna/lpa-sample/reference/kiv2_6.fasta",
                  mode = "complete"
                 ):
    """
    The function generate meta-informations of the VCFv4.2
    Args:
        file_format: String, the file format, default as VCFv4.2
        file_date: String, any datetime, when using default "today"
                   it will automatically generate the date
        source: String, the source, default is "lpa-analysis"
        reference: String, the reference, the default is the reference
                   of Coassin pipeline
        mode: String, the setting can accept "complete" and "simplified",
              under "complete" mode it will generate variant level(VL),
              total coverage(TC), and TypeB(TB)
    Return:
        String, the complete meta information
    """

    meta_information = \
f'''##fileformat={file_format}
##fileDate={datetime.datetime.now().strftime("%Y%m%d") if file_date=="today" else file_date}
##source={source}
##reference={reference}
##phasing=partial
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
'''

    if mode == "complete":
        meta_information +=\
"""##FORMAT=<ID=VL,Number=1,Type=Float,Description="Variant Level">
##FORMAT=<ID=TC,Number=1,Type=Integer,Description="Total Coverage">
##FORMAT=<ID=TB,Number=1,Type=String,Description="TypeB">
"""
    return meta_information

def generate_header(sample_id_list):
    """
    the function generate the header of the table
    Args:
        sample_id_list: list[str], the list of sample_id
    Return:
        the complete header
    """
    return "\t".join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT'] + sample_id_list)+"\n"

def new_index_VCF(l, i, mode):
    """
    a exception solver for list with exceed index, when exceed, return a placeholder
    Args:
        l: list[str]: a list of genotype literals
        i: int, the index
        mode: the mode in "complete" or 'simplified',
              determined the placeholder's format
    Return:
        the output, if the l[i] is invalid, return a placeholder
    """
    try:
        return l[i]
    except:
        if mode == "complete":
            return "0/0:::"
        else:
            return "0/0"

def new_index_LPA(l, i, mode, order):
    """
    a exception solver for list with invalid index, return the corresponding
    output based on the genotype literal in LPA mode
    Args:
        l: list[str]: a list of genotype literals
        i: int, the index
        mode: the mode in "complete" or 'simplified',
              determined the placeholder's format
    Return:
        the output, if the l[i] is invalid, return a placeholder
    """
    # if this item fit in order[i]'s pattern, i.e. order[i][0] 
    try:
        return [item[1] for item in l if item[0] == order[i][0]][0]
    except:
        if mode == "complete":
            return "0/0:::"
        else:
            return "0/0"

def comp_list(l):
    """
    comparison key between list with at most 2 element
    length is dominant, then compare first term, then second term
    helper function of assign_order as a comparison key
    """
    if len(l) == 1:
        return 100 + 10*l[0]
    else:
        return 200 + 10*l[0] + l[1]
    
def assign_order(l):
    """
    given a list of list, assign them a order given by comp_list,
    just because list of list cannot be hashed and sorted...
    """
    l = [eval(j) for j in set([str(i) for i in l])]
    l.sort(key = comp_list)
    return [[j,i] for i,j in enumerate(l)]


def generate_line(pos, REF, ALT, sample_literal, mode, format):
    """
    generating actual lines for each position <pos> with reference <REF> and
    alternatives <ALT> from <sample_literal>,
    Args:
        pos: int, the position
        REF, String, the reference on this position
        ALT, String, a comma-joined list of alternatives at this position
        sample_literal, list[list], all the genotype literals at this position
        mode, String, "complete" or "simplified"
        format, String, "VCF" or "LPA"
    Return:
        the actual lines in VCF file
    """
    CHROM = "6"
    POS = str(pos)
    ID = '.'
    QUAL = '.'
    FILTER = '.'
    if mode == "complete":
        FORMAT = "GT:VL:TC:TB"
    else:
        FORMAT = "GT"
    # in VCF mode, each ALT list all the alternatives
    # the genotype will use 1,2,3... in the literal and list them "from top to bottom"
    if format == "VCF":
        # get the max repeating line
        repeat = max([len(l) for l in sample_literal])
        if repeat == 1:
            # for position with only one line, get the first one in each sublist
            return "\t".join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, FORMAT]+[i[0] for i in sample_literal])+"\n"
        else:
            # for position multiple, list all possible results "top to bottom"
            # if a sample is running out of possibilities use placeholder
            return "\n".join(
                ["\t".join(
                    [CHROM, POS, ID, REF, ALT, QUAL, FILTER, FORMAT]+
                    [new_index_VCF(l,i,mode) for l in sample_literal])
                 for i in range(repeat)
                ]
            )+"\n"
    # in LPA mode, each alternatives will take a line
    # the genotype will use 1 in the literal and list them on the corresponding line
    else:
        # flatten the sample literal by 2 level then pick [0], 
        # i.e. get all the appeared pattern, assign an order to this list of list
        order = assign_order([l2[0] for l1 in sample_literal for l2 in l1 if len(l2) == 2])
        # get the total number of lines
        repeat = len(order)
        # each alternatives will be connected by commna
        ALT = ALT.split(",")
        ALT = [",".join([ALT[j-1] for j in i[0]]) for i in order]
        return "\n".join(
            ["\t".join(
                [CHROM, POS, ID, REF, ALT[i], QUAL, FILTER, FORMAT]+
                formatting_lines(sample_literal,i,mode,order))
             for i in range(repeat)
            ]
        )+"\n"


def formatting_lines(sample_literal,i,mode,order):
    #convert the lines into string
    body = str([new_index_LPA(l,i,mode,order) for l in sample_literal])
    # find all non 0 numbers, 
    numbers = list(set([i for j in re.findall(r'([0-9])/([0-9])', body) for i in j]))
    try:
        numbers.remove("0")
    except ValueError:
        pass
    if len(numbers) > 1:
        # it can be sorted although they are strings
        alt_list = sorted(numbers)
        body = body.replace(alt_list[0],"1")
        body = body.replace(alt_list[1],"2")
    elif len(numbers) == 1 and numbers[0] != "1":
        # it can be sorted although they are strings
        body = body.replace(numbers[0],"1")
    return eval(body)
        
    
def get_positional_data(df):
    """
    Extract the positional-specific data, df
    Args:
        df: pandas.DataFrame, which is expected to be the return of
        pandas.DataFrameGroupBy.get_group(), grouped by the position("Pos")
    Return:
        REF: String, the reference
        ALT: String, the comma-joined alternatives
        TC: 'Coverage-Total' column
        alt_list: list[string], the list format of ALT
    """
    # reference is unique for each position
    REF = df.Ref.unique()[0]
    # This list keeps the order of alternatives
    alt_list = list(set(list(df.Variant.unique()) +
               [alt for i in list(df["Major/Minor"].unique()) for alt in i.split("/")]))
    alt_list = [i for i in alt_list if i != REF]
    # the output format of alternatives
    ALT = ','.join(alt_list)
    # Total Coverage is unique for each position
    TC = df['Coverage-Total'].unique()[0]
    return REF, ALT, TC, alt_list

def formatting_alt(l):
    if (len(l) == 2 and l[1] == l[0]):
        return [l[0]]
    else:
        return sorted(l)
    
def _get_genotype_literal(row, ref, alt_list, format):
    '''
    generate the genotype information literal
    Args:
        row: pandas dataframe individual line, a variant record
        ref: str, the reference position
        alt_list: all the alternatives possible
        format: String, "VCF" or "LPA"
    Return:
        format == "VCF"
        String, a genotype output follows VCFv4.2 1.4.2 genotype field
                requirement, given the order from alt_list
        format == "LPA"
        list[list, str] a genotype output follows VCFv4.2 1.4.2 genotype field
                        requirement given the order from alt_list,
                        the list is the 1-based position in alt_list, 
                        if the length is 2, it means we have more than 1 alt in ref/alt
                        str is the genotype output
    '''
    gt_list = [ref]+alt_list
    gt_dict = {v: k for k, v in enumerate(gt_list)}
    major, minor = row['Major/Minor'].split('/')
    if format == "VCF":
        return f"{gt_dict[major]}/{gt_dict[minor]}"
    else:
        alt = [gt_dict[allele] for allele in [major, minor] if gt_dict[allele] != 0]
        return [formatting_alt(alt),f"{gt_dict[major]}/{gt_dict[minor]}"]

def _get_genotype(df, ref, alt_list, format):
    """
    a wrapper apply _genotype to each individual line of df in VCF format
    Args:
        df: pandas.DataFrame, a DataFrame of variant record
        ref: str, the reference position
        alt_list: all the alternatives possible
        format: "VCF" or "LPA"
    Return:
        dict{str:str} or dict{str, list[list[int], str]}
        genotype output follows VCFv4.2 1.4.2 genotype field
        requirement given the order from alt_list,
        first list[str] is the row index,
        second str is the genotype literal
    """
    # apply _get_genotype_literal to all rows
    return dict(
        df.apply(
            func = (
            lambda x: _get_genotype_literal(row = x,
                                            ref = ref,
                                            alt_list = alt_list,
                                            format = format
                                            )
                ),
            axis=1
            )
        )

def _extract_data(df, ref, alt_list, TC, mode, format):
    """
    generate proper data format from <df>
    Args:
        df: pandas.DataFrame object generated from a get_group method
            from a pandas.DataFrameGroupBy object, grouped by the position("Pos")
            then sample ID("SampleID")
        ref: String, the refernce at this position
        alt_list: list[str], the list of alternatives
        TC: String, the total coverage
        mode: "complete" or "simplified"
        format: "VCF" or "LPA"
    Return:
        format == "VCF"
        list[str], a genotype output follows VCFv4.2 1.4.2 genotype field
                requirement, given the order from alt_list
        format == "LPA"
        list[list[list[int], str]] a genotype output follows VCFv4.2 1.4.2 genotype field
                       requirement given the order from alt_list,
                       the int is the 1-based position in alt_list,
                       str is the genotype output

    """
    # extract the genotype notation string
    GT = _get_genotype(df = df, ref = ref, alt_list = alt_list, format = format)
    # when outputing VCF format
    if format == "VCF":
        if mode == "complete":
            # get the variant level information
            VL = dict(df['Variant-Level'])
            # get the typeB information
            TB = dict(df["TypeB"])
            # build the notation literals
            notation = [f"{GT[index]}:{VL[index]:.4f}:{TC}:{TB[index]}" for index in df.index]
        else:
            notation = [f"{GT[index]}" for index in df.index]
    # when outputing LPA format
    else:
        if mode == "complete":
            # get the variant level information
            VL = dict(df['Variant-Level'])
            # get the typeB information
            TB = dict(df["TypeB"])
            # build the notation literals
            notation = [[GT[index][0],f"{GT[index][1]}:{VL[index]:.4f}:{TC}:{TB[index]}"] for index in df.index]
        else:
            notation = [[GT[index][0],f"{GT[index][1]}"] for index in df.index]

    return notation

def get_sample_data(df_group_by_id, SampleID, ref, alt_list, TC, mode, format):
    """
    a wrapper for _extract_data, which generate a
    placeholder literal for the keyError
    (this SampleID doesn't have mutation at this position)
    Args:
        df: pandas.DataFrameGroupBy object, generated from a groupby method
          grouped by the position("Pos") then sample ID("SampleID")
        SampleID: String, the sample id
        ref: String, the refernce at this position
        alt_list: list[str], the list of alternatives
        TC: String, the total coverage
        mode: "complete" or "simplified"
        format: "VCF" or "LPA"
    Return:
        list[str]: a list with genotype literals inside
    """
    try:
        df = df_group_by_id.get_group(SampleID)
    except KeyError:
        if mode == "complete":
            return ["0/0:::"]
        else:
            return ['0/0']
    # if there is any duplicate use comma for now
    return _extract_data(
                        df = df,
                        ref = ref,
                        alt_list = alt_list,
                        TC = TC,
                        mode = mode,
                        format = format
                        )

def write_vcf(input_path = "data", bam_list = "*", output_path = "output.vcf", mode = "simplified", format = "LPA"):
    """
    Wrapper function for a complete process taking Coassin pipeline output as input
    and output VCF file at output_path

    Args:
        input_path: str, the path of input folder with all the Coassin Pipeline output inside the folder
        bam_list: list[str], if specified it will select the file in input_path
        output_path: str, the path of output path
        mode: str, in "complete" mode all the features will be printed,
              in "simplified" mode only the genotype will be generated
        format: str, in "VCF" mode it will list different genotypes together
                 in "LPA" mode it will list each individual genotypes on seperate lines
    output:
        a file at <output_path> following VCFv4.2 requirement
    """
    assert mode in ["complete", "simplified"]
    assert format in ['VCF', 'LPA']

    if bam_list == "*":
        bam_list = [dir for dir in next(os.walk(input_path))[1] if '.BQSR.recaled.bam' in dir]

    # Load the complete dataset into pandas
    df = data_load_wrapper(input_path, bam_list)

    # Write the output path
    write(path = output_path, content = generate_meta(mode = mode))

    # get all the SampleID(this ID is str, but it's still sorted numerically in python)
    sample_id_sorted_list = list([i for i in df.SampleID.sort_values().drop_duplicates()])

    # append the header
    append(path = output_path,
           content = generate_header(sample_id_list = sample_id_sorted_list))

    # get all the positions
    pos_sorted_list = list([i for i in df.Pos.sort_values().drop_duplicates()])

    # group the results by Positions("Pos"), sort = False to improve efficiency
    df_group_by_pos = df.groupby(['Pos'], sort = False)

    # iterate over the positions, to save the memory the result will be output line by line, i.e. by position
    for pos in pos_sorted_list:

        # get all the records on position <pos>
        df_group_pos = df_group_by_pos.get_group(pos)

        # get the positional data, which will not be affected or can be generated by individual Samples
        REF, ALT, TC, alt_list = get_positional_data(df_group_pos)

        # group the records by sample id ("SampleID"), sort = False to improve efficiency
        df_group_by_id = df_group_pos.groupby("SampleID", sort = False)

        # get data from each individual sample
        sample_literal = [get_sample_data(
            df_group_by_id = df_group_by_id,
            SampleID = sample_id,
            ref = REF,
            alt_list = alt_list,
            TC = TC,
            mode = mode,
            format = format
        ) for sample_id in sample_id_sorted_list]

        # complete the row with sample data
        append(path = output_path,
               content = generate_line(pos = pos,
                                       REF = REF,
                                       ALT = ALT,
                                       sample_literal = sample_literal,
                                       mode = mode,
                                       format = format))
