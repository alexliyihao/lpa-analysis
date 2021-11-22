import os
import pandas as pd
from . import settings
def add_prefix(parent_list, hispanic_list, i):
    """
    helper function for generate_bam_list, add correct "hispanic" path if necessary
    """
    if i in parent_list:
        return i
    elif i in hispanic_list:
        return os.path.join("hispanic",i)
    
def generate_bam_list(parent_path, csv_path):
    """
    the function generate all the relative bam files from a pandas.DataFrame
    Args:
        parent_path: str, the path taking all the coassin output
        csv_path: str, the dataframe to be analyzed
    Return:
        list[str]: the list of all relative bam files to parent path
    """
    df = pd.read_csv(csv_path, index_col = 0)
    parent_list = next(os.walk(parent_path))[1]
    hispanic_list = next(os.walk(os.path.join(parent_path, "hispanic")))[1]
    return [add_prefix(parent_list = parent_list, 
                       hispanic_list = hispanic_list, 
                       i = WES_ID+".BQSR.recaled.bam") 
            for WES_ID in df.WES_ID]

def generate_paths(ethnicity):
    """
    utility function generate all the related paths
    Args:
        ethnicity: str or int, the ethnicity node
    """
    assert ethnicity in [1,2,3,4, 'complete']
    INPUT_PATH = settings.COASSIN_OUTPUT
    BAM_LIST = generate_bam_list(parent_path = settings.COASSIN_OUTPUT, csv_path = settings.SNPS[ethnicity])
    VCF_NAME = f"ethnicity_{ethnicity}.vcf"
    VCF_PATH = os.path.join(settings.VCFs, VCF_NAME)
    CSV_NAME = f"ethnicity_{ethnicity}.csv"
    OUTPUT_PATH = os.path.join(settings.CSVs, CSV_NAME)
    return INPUT_PATH, BAM_LIST, VCF_NAME, VCF_PATH, CSV_NAME, OUTPUT_PATH