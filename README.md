# lpa-analysis

Code repo for *[Ancestry specific distribution of LPA Kringle IV-Type-2 genetic variants highlight associations to apo(a) copy number, glucose and hypertension](TODO: add links)*

## Environment

 - Coassin_pipeline's dependency is packaged in Coassin_pipeline/app_tar, for details, see the Readme file in Coassin_pipeline.
 - The project's python dependencies and releases are managed by [Poetry](https://python-poetry.org/). If you only works on codes, related information can be found at pyproject.toml under section tool.poetry.dependencies, For running details, see the Readme file in src/lpa_pipeline.
 - Some visualization code are written by R. The sessionInfo can be found at visualization/sessionInfo.txt

## Usage

- **Coassin_pipeline** included shell script and all the related example file used running pipeline from Dr. Coassin et.al's paper [**A comprehensive map of single base polymorphisms in the hypervariable LPA Kringle-IV-2 copy number variation region**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314250/), the [pipeline](https://github.com/genepi/lpa-pipeline) is running on [Cloudgene](http://www.cloudgene.io/) on a SGE or SLURM-managed cluster

- **src/lpa_pipeline** is the data analysis pipeline used in Python.

- **METAL** included the backup shell script running [METAL](https://github.com/statgen/METAL) 11.3.25 for meta analysis. These scripts are not required in the running

- For the VNTR estimating pipeline, please see the repo [vntrwrap](https://github.com/alexliyihao/vntrwrap). Those pipeline are list separately for they are also working for other projects.

## Documentation

- Most of the Python codes are documented in their source code
- A Sphinx-based autodoc is host [here](TODO: add the link)
- A copy of the autodoc can be found in docs/_build/html/index.html, which is saved on purpose.

## License

The repository was created by Yihao Li. It is licensed under the terms of the MIT license.

## Contact

If need help or can help with anything, please feel free to email Yihao Li (yl4326[at]cumc.columbia.edu) and CC both Dr. Gissette Reyes-Soffer (gr2104[at]cumc.columbia.edu) and Dr. Badri Vardarajan (bnv2103[at]cumc.columbia.edu)
