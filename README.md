# lpa-analysis

Code repo for research project with Dr. Gissette Soffer and Dr. Badri Vardarajan at Columbia Medical Center

## Environment

 - Coassin_pipeline's dependency is packaged in Coassin_pipeline/app_tar, for details, see the Readme file in Coassin_pipeline.
 - The project's dependencies and releases are managed by [Poetry](https://python-poetry.org/). If you only works on codes, related information can be found at pyproject.toml under section tool.poetry.dependencies, For running details, see the Readme file in src/lpa_pipeline.


## Usage

- **Coassin_pipeline** included shell script and all the related example file used running pipeline from Dr. Coassin et.al's paper [**A comprehensive map of single base polymorphisms in the hypervariable LPA Kringle-IV-2 copy number variation region**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314250/), the [pipeline](https://github.com/genepi/lpa-pipeline) is running on [Cloudgene](http://www.cloudgene.io/) on a SGE-managed cluster

- **src/pipeline** is the data analysis pipeline used in Python.

- **METAL** included the backup shell script running [METAL](https://github.com/statgen/METAL) 11.3.25 for meta analysis. These scripts are not required ion running

## Documentation

- A Sphinx-based api copy can be found in docs/_build/html/index.html
- TODO: host a documentation online(readthedocs.com, maybe?)

## License

The repository was created by Yihao Li. It is licensed under the terms of the MIT license.
