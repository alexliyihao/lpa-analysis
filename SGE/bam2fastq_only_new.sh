#########################################################################
## Section 1: Header                                                    #
#########################################################################

#!/bin/bash

# Specify name to be used to identify this run
#$ -N bam_to_fastq

# Email address (change to yours)
#$ -M example@columbia.edu

# Specify mailing options: b=beginning, e=end, s=suspended, n=never, a=aborted or suspended
#$ -m besa

# This sets the task range in the array from 1 to 5 with a step size of 1
#$ -t 1-5:1

# Change directory to the current
#$ -cwd

# Specify that bash shell should be used to process this script
#$ -S /bin/bash

# Specify the output file
#$ -o ~/bam2fastq_run/$JOB_NAME_$TASK_ID.out

# Specify the error file
#$ -e ~/bam2fastq_run/$JOB_NAME_$TASK_ID.err



#########################################################################
## Section 2: Transfering Data from current directory to /tempdata      #
#########################################################################

# Specify the filenames
INPUTFILES=(washei71450.BQSR.recaled.bam washei71453.BQSR.recaled.bam washei71451.BQSR.recaled.bam washei71454.BQSR.recaled.bam washei71452.BQSR.recaled.bam)

INPUTPATH="[input_path]"

# Name the Project directory
PROJECT_SUBDIR="ArrayJob_bam2fastq"

# Pull data file name from list defined above according to job id
INPUTFILENAME="$INPUTPATH/${INPUTFILES[$SGE_TASK_ID - 1]}"

# Specify working directory under the /tempdata
WORKING_DIR="$USER/tempdata/$PROJECT_SUBDIR-$SGE_TASK_ID"

# If working directory does not exist, create it
# The -p means "create parent directories as needed"
if [ ! -d "$WORKING_DIR" ]; then
    mkdir -p $WORKING_DIR
fi

# Specify destination directory (this will be subdirectory of your user directory in the archive)
DESTINATION_DIR="$USER/$PROJECT_SUBDIR/$PROJECT_SUBDIR-$SGE_TASK_ID"

# If destination directory does not exist, create it
# The -p in mkdir means "create parent directories as needed"
if [ ! -d "$DESTINATION_DIR" ]; then
    mkdir -p $DESTINATION_DIR
fi

# Copy the input data to the working directory
cp $INPUTFILENAME $WORKING_DIR/${INPUTFILES[$SGE_TASK_ID - 1]}

# Copy the app to the working directory
cp ~/bam2fastq.tgz $WORKING_DIR/bam2fastq.tgz

#########################################################################
## Section 3:Executing the program                                      #
#########################################################################

# navigate to the working directory
cd $WORKING_DIR

# Prepare a folder for the output
mkdir fastq

# deploy the app
tar -xvf bam2fastq.tgz
cd bam2fastq
make

# Run the program
./bam2fastq -o ../fastq/${INPUTFILES[$SGE_TASK_ID - 1]}_output#.fastq ../${INPUTFILES[$SGE_TASK_ID - 1]}


#########################################################################
## Section 4: Copy the result to destination                            #
#########################################################################
cd ../fastq
cp * ~/$DESTINATION_DIR

# clear the working directory
cd ../../../..
rm -rf $WORKING_DIR
