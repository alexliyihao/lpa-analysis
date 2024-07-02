# Coassin Pipeline

Shell script and all the related example file used running pipeline from Dr. Coassin et.al's paper [**A comprehensive map of single base polymorphisms in the hypervariable LPA Kringle-IV-2 copy number variation region**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314250/), the [pipeline](https://github.com/genepi/lpa-pipeline) is running on [Cloudgene](http://www.cloudgene.io/) on a SGE or SLURM-managed cluster

## Environment
Due to a scheduler transition in CUIMC Neurology HPC cluster, both scripts for Sun Grid Engine(SGE) and SLURM are provided. Which is running on a Debian environment, running the following pipelines:

## Steps
0. Input are .bam files from Washington Heights-Hamilton Heights-Inwood Columbia Aging Project (WHICAP) dataset
1. bam2fastq https://github.com/jts/bam2fastq
2. bwa-mem https://github.com/lh3/bwa
3. the pipeline from Dr. Coassin et.al's paper **A comprehensive map of single base polymorphisms in the hypervariable LPA Kringle-IV-2 copy number variation region** (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314250/), the pipeline can be find at https://github.com/genepi/lpa-pipeline running on cloudgene http://www.cloudgene.io/
4. Outputing tab-delimiting .txt files

## Note:
 - Please notice that although all the intermediate results are deleted once finished its task, this pipeline is still extremely disk-storage-demanding when using complete exome data(the fastq file will be ~120GB per job at its peak, although the final output is in KBs), and either 1. prepare enough storage spaces for a smooth running. 2. consider using *REGION=""* option (see line 45-46 and the comments above).

 - For some reason, the pipeline install and delete bam2fastq, bwa and cloudgene in each individual sub-task, it's completely okay to simply call them from somewhere else by modifying section 3 **except** the cloudgene (each job extract cloudgene's latest job finished, I cannot figure out a easy way to extract corresponding job for each bam file).

 - In the original pipeline the bam file have a sub-directory (see data_inflow/bam_list.txt), you may want to delete some lines in the .sh script.
