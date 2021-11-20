The script translating Coassin pipeline's output to VCF format
```
usage: lpa_vcf_translate.py [-h] [-i INPUT_PATH] [-o OUTPUT_PATH]
                            [-bl BAM_LIST] [-f FORMAT] [-m MODE]

optional arguments:

  -h, --help            show this help message and exit
  -i INPUT_PATH, --input_path INPUT_PATH
                        The path of folder storing all the output folders
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output file
  -bl BAM_LIST, --bam_list BAM_LIST
                        The list of bam files to be included, in python list
                        format, * will search for all compatible folders
  -f FORMAT, --format FORMAT
                        The format of output, LPA or VCF
  -m MODE, --mode MODE  The output mode, complete or simplified
```
