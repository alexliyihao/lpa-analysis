# LPA pipeline

The code for statistical analysis

## Import

The pipeline is written as Python modules, a demo JupyterLab IPython Notebook(.ipynb) is provided in `docs/_build/html/index.html`, due to some version issue we cannot host it on readthedoc.com.

The tqdm are included in the code, in case there's anything wrong with tqdm displaying, our running environment is using jupyterlab==3.2.4 with widgetsnbextension==3.6.4.

- As python module: the usage is provided in module docstring.
  ```
  from lpa_pipeline import <filename>
  print(<filename>.__doc__)
  ```
- (for some of the scripts) As command line script: the usage is provided by argparse description. But not tested
  ```
  python <filename>.py -h
  ```

## Running

1. Extract the information from Coassin's Pipeline
    This step takes the output of SGE task in `Coassin_pipeline/coassin_pipeline.sh` as input. Pipeline Code running the following steps:

    encodings.py:
     - encode the Coassin's output into one-hot(0/1) encoded carrier information.
     - You may want to modify _extract_ID() method for your data, it takes in the string name of output folder and return a new string as ID.

    locus_collector.py
     - extract the mylocus info from Coassin's output

2. Filter the SNPs encoded
This step takes the output of encodings.py as input. Pipeline Code running the following steps:

    snps_filter.py
    - filter the SNPs with bad reading qualities in the dataset

3. Association analysis

    This step takes the following as inputs:
   1. the output of step 2, i.e. filtered SNPs table
   2. Other exogenous variables table
   3. target variable table
   All of which should be row-aligned

    Association.py
     - Iteratively do the following:
     ```
    for y in target variable table:
        for snp in filtered SNPs table:
          (if extra_iterate_on specified) for each ethnicity:
            run regression "y ~ snp + Other exogenous variables table" (R format here)
    ```
     - For its iterative nature and various pipeline setting, this one only runs as Python module.

4. METAL meta-analysis

    This step takes the output of step 3 as inputs

    METAL_toolkit.py
    - copy the result of Association.py to a new folder
      - The purpose for this copying is to keep an untouched version of data.
    - Run METAL internally, all the scripts are generated by the code.
    - Aggregate all the results
    - The output is a .xlsx file rather than a csv file - for it will usually be used for exploratory operations

5. Post Processing

    Takes the output os step 4 as inputs

    post_processing.py
    - Compute the FDR-adjusted p-value
    - Correct the ``METAL``'s direction discrepancy
    - Appended the locus information necessary

## Visualizations

1. freq_table_generator.py
   - Compute the SNPs frequency in each class defined in <class_variable>, used for Supplemental table S1

2. Appearance table generator.py
   - A generator computing appearance absolute frequency of SNP carrier by race and ethnicity, used for Supplemental Figure S2

3. haploview_generator.py
   - A class generating output for Haploview4.2 linkage format, used for [Haploview](https://www.broadinstitute.org/haploview/haploview) and supplemental table S3D
