## Section 1: Header ----------------------------------------------------

#!/bin/bash

# Specify name to be used to identify this run
#$ -N lpa-analysis-LASSO

# Email address (change to yours)
#$ -M yl4326@columbia.edu

# Specify mailing options: b=beginning, e=end, s=suspended, n=never, a=aborted or suspended
#$ -m besa

# This sets the task range in the array
#$ -t 1-132:1

# Change directory to the current
#$ -cwd

# Specify that bash shell should be used to process this script
#$ -S /bin/bash

# Specify the output file
#$ -o /mnt/mfs/hgrcgrid/homes/$USER/logs/$JOB_NAME_$TASK_ID.outerr

# Specify the error file
#$ -e /mnt/mfs/hgrcgrid/homes/$USER/logs/$JOB_NAME_$TASK_ID.outerr




## Section 2: Transfering Data & Path Settings ----------------------------------------------------
echo "Setting: Current Working directory is $PWD, current task id is $SGE_TASK_ID"

# Specify the data inflow and home directory
HOMEPATH="/mnt/mfs/hgrcgrid/homes/$USER"

# Specifying the name of the data
cd $HOMEPATH
module load CONDA

# Name the project subdirectory
PROJECT_SUBDIR="lasso_Nov27"

# Specify working directory under the /tempdata
WORKING_DIR="$PROJECT_SUBDIR/tempdata/$PROJECT_SUBDIR-$SGE_TASK_ID"
echo "Setting: The temporary working directory is $WORKING_DIR"

# If working directory does not exist, create it
# The -p means "create parent directories as needed"
if [ ! -d "$WORKING_DIR" ]; then
    mkdir -p $WORKING_DIR
fi



## Section 3: Deploying the Program ----------------------------------------------------
# Copy the input data to the working directory
cp $HOMEPATH/LASSO.tar $HOMEPATH/$WORKING_DIR/LASSO.tar
echo "Deployment: copied LASSO.tar from $HOMEPATH to $HOMEPATH/$WORKING_DIR/LASSO.tar"

# Navigate to the working directory
cd $HOMEPATH/$WORKING_DIR

# Deploy the bam2fastq
# Extract the file in slient mode
tar -xf LASSO.tar
echo "Deployment: Extracted LASSO.tar"
# navigate back to the working directory
cd LASSO

## Section 4:Executing the Program ----------------------------------------------------
echo "Run: Attempting to run Lasso"
conda run -n lasso python Lasso.py -i `expr $SGE_TASK_ID - 1`
echo "Run: Finished the Lasso"
cd $HOMEPATH
rm -rf $WORKING_DIR
