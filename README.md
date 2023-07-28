# lpa-analysis

Code repo for research project with Dr. Gissette Soffer and Dr. Badri Vardarajan at Columbia Medical Center

## Usage

- **Coassin_pipeline** included shell script and all the related example file used running pipeline from Dr. Coassin et.al's paper [**A comprehensive map of single base polymorphisms in the hypervariable LPA Kringle-IV-2 copy number variation region**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314250/), the [pipeline](https://github.com/genepi/lpa-pipeline) is running on [Cloudgene](http://www.cloudgene.io/) on a SGE-managed cluster

- **environment** provides [Anaconda](https://www.anaconda.com/) environment description files for fast virtual environment building.

- **pipeline** is the data analysis pipeline used in Python.
```{include} src/lpa_pipeline/README.md
```


- **METAL** included the backup shell script running [METAL](https://github.com/statgen/METAL) 11.3.25 for meta analysis. These scripts are not required ion running

## License

The repository was created by Yihao Li. It is licensed under the terms of the MIT license.
