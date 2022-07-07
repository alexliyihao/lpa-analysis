## Section 1: Header ----------------------------------------------------

#!/bin/bash

# Specify name to be used to identify this run
#$ -N lpa-analysis

# Email address (change to yours)
#$ -M yl4326@columbia.edu

# Specify mailing options: b=beginning, e=end, s=suspended, n=never, a=aborted or suspended
#$ -m besa

# This sets the task range in the array
#$ -t 1-3917:1

# This set concurrent job limit(limit to 4 for the disk usage)
# -tc 4(not applied here)

# Change directory to the current
#$ -cwd

# Specify that bash shell should be used to process this script
#$ -S /bin/bash

# Specify the output file
#$ -o /mnt/mfs/hgrcgrid/shared/LPA_analysis/logs/$JOB_NAME_$TASK_ID.outerr
# If you prefer splitting out and error
# -o /mnt/mfs/hgrcgrid/shared/LPA_analysis/logs/$JOB_NAME_$TASK_ID.out

# Specify the error file
#$ -e /mnt/mfs/hgrcgrid/shared/LPA_analysis/logs/$JOB_NAME_$TASK_ID.outerr
# If perfer splitting
# -e /mnt/mfs/hgrcgrid/shared/LPA_analysis/logs/$JOB_NAME_$TASK_ID.err




## Section 2: Transfering Data & Path Settings ----------------------------------------------------
echo "Setting: Current Working directory is $PWD, current task id is $SGE_TASK_ID"

# Specify the data inflow and home directory
INPUTPATH="/mnt/mfs/hgrcgrid/data/whicap/WHICAP_WES/BAM/washeiDragenBamsList/washeiBamsUpdate2/BQSR/bqsrRGbam"
HOMEPATH="/mnt/mfs/hgrcgrid/shared/LPA_analysis"

# Specifying the name of the data
cd $HOMEPATH
readarray -t INPUTFILES < data_inflow/bam_list.txt 
echo "Setting: ${#INPUTFILES[@]} files detected"
echo "Setting: working on ${INPUTFILES[$SGE_TASK_ID - 1]}"

# Pull data file name from list defined above according to job id
INPUTFILENAME="$INPUTPATH/${INPUTFILES[$SGE_TASK_ID - 1]}"

# Deal with the additional /hispanic, preventing additional folders created
FILENAME_CLEANED="`echo "${INPUTFILES[$SGE_TASK_ID - 1]}" | sed 's/hispanic\///g'`"
echo "Setting: For input $INPUTFILENAME, cleaned into $FILENAME_CLEANED"

# Name the project subdirectory
PROJECT_SUBDIR="pipeline_output"

# Specify working directory under the /tempdata
WORKING_DIR="$PROJECT_SUBDIR/tempdata/$PROJECT_SUBDIR-$SGE_TASK_ID"
echo "Setting: The temporary working directory is $WORKING_DIR"

# If working directory does not exist, create it
# The -p means "create parent directories as needed"
if [ ! -d "$WORKING_DIR" ]; then
    mkdir -p $WORKING_DIR
fi

# Specify destination directory, if the original file is in /hispanic, it should still be in /hispanic
DESTINATION_DIR="$PROJECT_SUBDIR/${INPUTFILES[$SGE_TASK_ID - 1]}"
echo "Setting: The destination directory is $DESTINATION_DIR"

# If destination directory does not exist, create it
# The -p in mkdir means "create parent directories as needed"
if [ ! -d "$DESTINATION_DIR" ]; then
    mkdir -p $DESTINATION_DIR
fi




## Section 3: Deploying the Program ----------------------------------------------------
# Copy the input data to the working directory
cp $INPUTFILENAME $HOMEPATH/$WORKING_DIR/$FILENAME_CLEANED
echo "Deployment: copying original data $INPUTFILENAME to $HOMEPATH/$WORKING_DIR/$FILENAME_CLEANED"

# Copy the materials to the working directory
# Bam2fastq
cp $HOMEPATH/app_tar/bam2fastq.tgz $HOMEPATH/$WORKING_DIR/bam2fastq.tgz
# bwa
cp $HOMEPATH/app_tar/bwa.tgz $HOMEPATH/$WORKING_DIR/bwa.tgz
# cloudgene
cp $HOMEPATH/app_tar/cloudgene.tgz $HOMEPATH/$WORKING_DIR/cloudgene.tgz
# Coassin's reference and annotations
cp $HOMEPATH/app_tar/coassin_pipeline_data.tgz $HOMEPATH/$WORKING_DIR/coassin_pipeline_data.tgz

# Navigate to the working directory
cd $HOMEPATH/$WORKING_DIR

# Deploy the bam2fastq
# Extract the file in slient mode
tar -xf bam2fastq.tgz
# Navigate into bam2fastq folder
cd bam2fastq
# Make bam2fastq in slient mode, inhibiting all warnings
# (although there are still some in the outerr file)
make -ws
echo "Deployment: bam2fastq deployed"

# navigate back to the working directory
cd $HOMEPATH/$WORKING_DIR

# delete the tar of bam2fastq
rm bam2fastq.tgz
echo "Deployment: bam2fastq installation package cleaned"

# Deploy the bwa-mem
# Extract the file in silent mode
tar -xf bwa.tgz
# Navigate into the bwa folder
cd bwa
# Make bwa in silent mode
make -ws
echo "Deployment: bwa deployed"

# Navigate back to the working directory
cd $HOMEPATH/$WORKING_DIR

# Delete the tar of bwa
rm bwa.tgz
echo "Deployment: bwa installation package cleaned"

# Deploy the cloudgene with Coassin Pipeline
tar -xf cloudgene.tgz
echo "Deployment: cloudgene deployed"

# Delete the tar of cloudgene
rm cloudgene.tgz
echo "Deployment: cloudgene installation package cleaned"

# Extract the necessary data
tar -xf coassin_pipeline_data.tgz
echo "Deployment: necessary data installed"

# Eelete the tar of necessary data
rm coassin_pipeline_data.tgz
echo "Deployment: necessary data installation package cleaned"




## Section 4:Executing the Program ----------------------------------------------------

# Navigate to the working directory
cd $HOMEPATH/$WORKING_DIR

# Prepare a folder for the fastq output
mkdir fastqs

# Run the program
cd bam2fastq
# Pass the .bam file, generate a <bam_name>_output_1.fastq and <bam_name>_output_2.fastq in fastqs folder
./bam2fastq -o $HOMEPATH/$WORKING_DIR/fastqs/${FILENAME_CLEANED}_output#.fastq \
              $HOMEPATH/$WORKING_DIR/$FILENAME_CLEANED
echo "Execution: bam2fastq finished"

# Get back to the working directory
cd $HOMEPATH/$WORKING_DIR

# Remove the original file
rm $FILENAME_CLEANED

# Prepare a folder for the bam output
mkdir bams

# Run bwa
cd bwa
# bwa need to index the reference provided by Coassin first
./bwa index $HOMEPATH/$WORKING_DIR/coassin_pipeline_data/kiv2_6.fasta

# Align <bam_name>_output_1.fastq and <bam_name>_output_2.fastq against Coassin reference,
# -v 1 only print errors, 0-4 will increase the verbosity, default 3 will give 10K+ line output
./bwa mem -v 1 \
  $HOMEPATH/$WORKING_DIR/coassin_pipeline_data/kiv2_6.fasta \
  $HOMEPATH/$WORKING_DIR/fastqs/${FILENAME_CLEANED}_output_1.fastq \
  $HOMEPATH/$WORKING_DIR/fastqs/${FILENAME_CLEANED}_output_2.fastq \
  | gzip -3 > $HOMEPATH/$WORKING_DIR/bams/$FILENAME_CLEANED
echo "Execution: bwa-mem finished at $HOMEPATH/$WORKING_DIR/bams/$FILENAME_CLEANED"

# Navigate to the working directory
cd $HOMEPATH/$WORKING_DIR

# Delete the intermediate .fastq files
rm -rf fastqs
echo "Execution: fastq files cleaned"

# Navigate into the cloudgene
cd cloudgene
# Run Coassin pipeline with bam input in bams folder,
# refernce, annotate base, and annotate region are provided by Coassin
./cloudgene run lpa-mutation-server \
  --input $HOMEPATH/$WORKING_DIR/bams \
  --archive $HOMEPATH/$WORKING_DIR/coassin_pipeline_data/kiv2_6.fasta \
  --annotateBase $HOMEPATH/$WORKING_DIR/coassin_pipeline_data/typeb_annotation.txt \
  --annotateRegion $HOMEPATH/$WORKING_DIR/coassin_pipeline_data/maplocus_v3.txt
echo "Execution: cloudgene finished"

cd $HOMEPATH/$WORKING_DIR
# Delete the intermediate .bam files
rm -rf bams
echo "Execution: bam files cleaned"




## Section 5: Copy the Result to Destination ----------------------------------------------------

# Navigate to Cloudgene, where the result is saved
cd cloudgene

# Find the latest "job-[date]" files available
# Find a directory whose name in "job-<whatever>" format,
# Reversely sorted by time, pick the first one, then deal with additional spaces
recent_job=$(find -type d \
                  -name "job-*" \
                  -printf '%T+ %p\n' \
                  | sort -r \
                  | head -1 \
                  | cut -d ' ' -f2 -)
echo "Finishing: The most recent job $recent_job found in Cloudgene"

# Navigate to the job folder
cd $recent_job

# Copy all the outputs in job founded to DESTINATION_DIR
cp -r * $HOMEPATH/$DESTINATION_DIR
echo "Finishing: copying from $PWD to $DESTINATION_DIR"

# Clear the working directory
cd $HOMEPATH
rm -rf $WORKING_DIR
echo "Finishing: Complete, cleaning the working directory"
