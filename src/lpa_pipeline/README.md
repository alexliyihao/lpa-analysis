For iteration purpose, the pipeline is written as Python modules, a demo JupyterLab IPython Notebook(.ipynb) and it's pdf preview are provided in `docs` folder. A command line version is provided as well.

All the <filename>.py file can be used either except association.py:
- As python module: the usage is provided in module docstring.
  ```
  from lpa_pipeline import <filename>
  print(<filename>.__doc__)
  ```
- As command line script: the usage is provided by argparse description.
  ```
  python <filename>.py -h
  ```

1. Extract the information from Coassin's Pipeline
    This step takes the output of SGE task in `Coassin_pipeline/coassin_pipeline.sh` as input. Pipeline Code running the following steps:

    encodings.py:
     - encode the Coassin's output into one-hot(0/1) encoded carrier information.
     - You may want to modify _extract_ID() method for your data, it takes in the string name of output folder and return a new string as ID.

    locus_collector.py
     - extract the mylocus info from Coassin's output


2. Filter the SNPs encoded
    This step takes the output of encodings.py as input. Pipeline Code running the following steps:

    snps_filter.py - filter the SNPs with bad reading qualities in the dataset

3. Association analysis
    This step takes the following as inputs:
     1. the output of step 2, i.e. filtered SNPs table
     2. Other exogenous variables table
     3. target variable table
     All of which should be row-aligned

     Association.py
     - Iteratively do the following:
       for y in target variable table:
          for snp in filtered SNPs table
              (if extra_iterate_on specified, split the group by ethnicity)
              run regression "y ~ snp + Other exogenous variables table" (R format here)
     - For its iterative nature and various pipeline setting, this one only runs as Python module.
       If you are command line user, please create your own setting from demo and association.py,
       save them into a python script, and run python <your_script>.py
     You may want to check both
     ```
     from lpa_pipeline import association
     print(association.__doc__)
     ```
     and
     ```
     snp_asso = SNPAssociation()
     print(snp_asso.fit.__doc__)
     ```
4. METAL meta-analysis
    This step takes the output of step 3 as inputs

    METAL_toolkit.py
    - copy the result of Association.py to a new folder
      - The purpose for this copying is to keep a untouched version of data.
    - Run METAL internally, all the scripts are generated by the code.
    - Aggregate all the results
    - The output is a .xlsx file rather than a csv file - for it usually will be used for exploratory operations

5. Frequency table generator
    This step takes the following as inputs:
     1. the output of step 2, i.e. filtered SNPs table
     2. Other exogenous variables table, which has a column(<class_variable>) indicating a classification of the whole population
     These two tables should be row-aligned.

     freq_table_generator.py
       - Compute the SNPs frequency in each class defined in <class_variable>