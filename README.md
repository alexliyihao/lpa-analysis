# lpa-analysis
Code repo for research project with Dr. Gissette Soffer and Dr. Badri Vardarajan at Columbia Medical Center

- Coassin_pipeline included shell script and all the related example file used running pipeline from Dr. Coassin et.al's paper [**A comprehensive map of single base polymorphisms in the hypervariable LPA Kringle-IV-2 copy number variation region**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314250/), the [pipeline](https://github.com/genepi/lpa-pipeline) is running on [Cloudgene](http://www.cloudgene.io/) on a SGE managed cluster


- VCF is the python script used to convert lpa output to VCFv4.2 files

- data_analysis is the script and output used to analysis the association between genotypes and the SNPs

- EIGEN included the python and shell script running [EIGENSOFT](https://github.com/DReichLab/EIG) 6.1.4 for SmartPCA algorithm
  For the project is running on a cluster importing EIGEN with environment modules(module loads),
  the python and shell scripts' interaction may seems not quite clear, a readme file is provided

- METAL included the python and shell script running [METAL](https://github.com/statgen/METAL) 11.3.25 for meta analysis
  For the project is running on a cluster importing METAL with environment modules(module loads),
  the python and shell scripts' interaction may seems not quite clear, a readme file is provided
