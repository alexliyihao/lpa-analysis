import pandas as pd

def get_vcf_names(vcf_path):
    """
    the function get column names from a VCF file, helper function for read_vcf
    Args:
        vcf_path: the path of the files
    Return:
        vcf_names: list[str], the column name of the vcf file
        i: int, the number of lines where the headers are at.
    """
    i = 0
    with open(vcf_path, "rt") as file:
        # keep reading the lines until finding the proper line
        for line in file:
            # the line start with "#CHROM"
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
            i += 1
    return vcf_names, i

def read_vcf(path):
    """
    the function reading a simpfies form VCF files into useable csv file
    Args:
        path: str, the path to the file
    Return:
        pandas.DataFrame, the tranposed pandas dataframe can be used in the next step
    """
    # get the column name of the vcf file
    names, i = get_vcf_names(path)
    # read the csv file
    vcf = pd.read_csv(path, delim_whitespace=True, header=None, names=names)
    # drop the information lines, i+1 is to drop the header line
    vcf.drop(range(i+1),axis = 0, inplace = True)
    # the last column name always have a '\n', drop it
    vcf.rename(columns={vcf.columns[-1]:vcf.columns[-1][:-1]}, inplace = True)
    # position-ref/alt is the only keys we will use
    vcf['pos-ref/alt'] = vcf.POS+"-"+vcf.REF+"/"+vcf.ALT
    # Drop not-related columns
    vcf.drop(["#CHROM", "ID", "QUAL", "FILTER", "FORMAT", "REF", 'ALT', 'POS'], axis = 1, inplace = True)
    # use pos-ref/alt as the index
    vcf.set_index("pos-ref/alt", inplace = True)
    return vcf.transpose()

def encode(df):
    """
    The function one-hot encoded the DataFrame
    Args:
        df: pandas.DataFrame, the output from read_vcf
    Return:
        df: pandas.DataFrame, the DataFrame which is one-hot-encoded
    """
    df = pd.get_dummies(df, drop_first = True, prefix_sep = "-")
    # from alphanumerical order all the 0/0 should be dropped, just in case
    assert(len([i for i in df.columns if "0/0" in i]) == 0)
    return df
    