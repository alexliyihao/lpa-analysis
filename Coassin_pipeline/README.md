The script is for deployment on Grid Engine(SGE) running on a Debian environment, running the following pipelines:

0. Input are .bam files from Washington Heights-Hamilton Heights-Inwood Columbia Aging Project (WHICAP) dataset
1. bam2fastq https://github.com/jts/bam2fastq
2. bwa-mem https://github.com/lh3/bwa
3. the pipeline from Dr. Coassin et.al's paper **A comprehensive map of single base polymorphisms in the hypervariable LPA Kringle-IV-2 copy number variation region** (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314250/), the pipeline can be find at https://github.com/genepi/lpa-pipeline running on cloudgene http://www.cloudgene.io/
4. Outputing tab-delimiting .txt files

Please notice that although all the intermediate results are deleted once finished its task, this pipeline is still extremely disk-storage-consuming (at most 120GB per job at its peak, although the final output is in KBs), and please prepare enough storage spaces for a smooth running.
