#!/bin/bash

source globals.sh

# This is the app being submitted.
APPNAME="tomanaRP4k"

# This is the prefix given to every output directory created by running one of 
# the all-trks jobs. This is defined by Luiz within "submit-condorRP.csh".
JOB_DIR_PREFIX="submit-"

# All of the inputs to "submit-condorRP.csh" are located here along with a 
# corresponding job_id.
JOB_CSV="all-trks.csv"

# This is the format of jobs.
JOB_FILE_FORMAT="$INFO_PREFIX job_id,job_name,inputfile"

if (( $# != 1 )); then
    echo "****************************************"
    echo "Run one of the jobs from all-trks."
    echo "The job_ids are numbered 1-14."
    echo ""
    echo "Usage:"
    echo "./runjob.sh <job_id>"
    echo "****************************************"
    exit 1
fi

JOB_ID=$1

echo -e $JOB_FILE_FORMAT

# Parse the job csv.
while IFS=, read -r id jobname inputfile
do
    echo -e "$INFO_PREFIX $id,$jobname,$inputfile"

    if [ "$id" = "$JOB_ID" ]; then

        # Set this to indicate that the job_id entered was actually found.
        # Used purely for the output seen at the end of this file.
        found_job=1

        # This is the directory created when running the job_id.
        JOB_DIR="$JOB_DIR_PREFIX$jobname"

        # Check if the job directory already exists, indicating that the job has
        # already been ran.
        if [ -d $JOB_DIR ]; then
            print_error "Directory $JOB_DIR already exists!"
            print_error "Make sure the job hasn't already ran."
            exit 1
        fi        

        # This is the output file created by this script.
        OUTPUT_FILE="output_$jobname.root"

        # Ensure the output file doesn't already exist as well.
        if [ -f "$OUTPUT_FILE" ]; then
            print_error "$OUTPUT_FILE exists already!"
            print_error "Have you already ran this job?"
            exit 1
        fi
        

        # submit the jobs.
        execute "csh ./submit-condorRP.csh $APPNAME $jobname $inputfile"

        exit 0
    fi
done < $JOB_CSV


# If this value is not set, tell the user the job_id was not found.
if [ -z ${found_job+x} ]; then echo "Unable to find job $JOB_ID."; fi
