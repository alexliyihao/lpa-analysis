The python and shell script running [EIGENSOFT](https://github.com/DReichLab/EIG) 6.1.4 for SmartPCA algorithm

The project is running on a cluster importing EIGENSOFT with Linux environment modules(module loads EIGEN/6.1.4), so it's a little bit inconvenient that the python instance cannot run the pipeline from start to the end.

1. EIGENToolkit:

    EIGENToolkit.py provides a EIGENToolkit class which takes several table as its input, and provides the following APIs:

    - Initialize:
        ```
        from EIGEN import EIGENToolkit
        etk = EIGENToolkit.EIGENToolkit(chromosome = 6,
                                        genetic_position = 0,
                                        physical_position_base = 161033785,
                                        verbose = 1)
        ```
        - chromosome, genetic_position, physical_position_base: check EIGENSOFT's documentation
        - verbose: 0 or 1, if 1, when new file is generated and written, a notice will be printed

    - One-line preprocessing APIs:

        - etk.to_eigen(snp, ethnicity_from_pca, gender, output_path, filename)

          Prepare a EIGENSTRAT format folder at output_path. All the files will be named as filename.ext
          The following three pd.DataFrame/Series should use individual sample as rows
          snp: pd.DataFrame, the snp integer-encoding table,
                       mutations as columns and individual sample as rows
          ethnicity_from_pca: pd.DataFrame, the ethnicity of each individual
          gender: pd.Series, the gender of each individuals

2. Running the EIGENSOFT SmartPCA:

    ```
    module loads EIGEN/6.1.4
    smartpca -p example.pca.par > example.log
    ```
    the script can be generated from *etk.to_eigen(path)* or an example *complete.pca.par* is provided as well. Please do not use the argument flag API
    ```
    smartpca.perl -i test.eigenstratgeno -a test.snp -b test.ind -k 3 -o test.pca -p test.plot -e test.eval -l test.log
    ```
    for it doesn't provide the PC weight functionality.


3. Reading helper:

    - etk.read_PC_weight(path): read a .snpweight output file and tidy the format

    - etk.read_PCA_result(path): read a .pca.evec output file and tidy the format
