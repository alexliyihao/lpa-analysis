The script is for deployment on Grid Engine(SGE) running on a Debian environment, calling the following pipelines:

- Input are .bam files from Washington Heights-Hamilton Heights-Inwood Columbia Aging Project (WHICAP) dataset

- bam2fastq https://github.com/jts/bam2fastq

- bwa-mem https://github.com/lh3/bwa

- the pipeline from Dr. Coassin et.al's paper **A comprehensive map of single base polymorphisms in the hypervariable LPA Kringle-IV-2 copy number variation region** (http://www.jlr.org/content/early/2018/11/09/jlr.M090381.full.pdf+html), the pipeline can be find at https://github.com/genepi/lpa-pipeline running on cloudgene http://www.cloudgene.io/

- Outputing tab-delimiting .txt files
