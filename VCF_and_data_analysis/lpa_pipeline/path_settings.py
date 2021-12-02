ROOT = "/mnt/mfs/hgrcgrid/shared/LPA_analysis"

DATAMETA = {
    "snps":{
        "path":"snps",
        "name_format": "pheno_snps_ethnicity_*.csv"
    },
    "genotypes":{
        "path":"genotypes",
        "name_format": "geno_ethnicity_*.csv"
    },
    "vcfs":{
        "path":"vcfs",
        "name_format": "ethnicity_*.vcf"
    }
}

COASSIN_OUTPUT = f"{ROOT}/coassin_pipeline/pipeline_output"

ANALYSIS_RESULT_PATH = "data_analysis_result"

CLASS_ALIAS = {"eu": 1,
               "af": 2,
               "hisp":3,
               "others":4,
               "complete": "complete"}

WHICAP_SOURCE = "/mnt/mfs/hgrcgrid/data/whicap/WHICAP_WES/BAM/washeiDragenBamsList/washeiBamsUpdate2/BQSR/bqsrRGbam"

