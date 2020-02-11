#!/bin/bash

# This is the prefix given to every output directory created by running one of 
# the all-trks jobs. This is defined by Luiz within "submit-condorRP.csh".
JOB_DIR_PREFIX="submit-"

# All of the inputs to "submit-condorRP.csh" are located here along with a 
# corresponding job_id.
JOB_CSV="all-trks.csv"

if (( $# != 1 )); then
    echo "****************************************"
    echo "Hadd one of the jobs from all-trks."
    echo "The job_ids are numbered 1-14."
    echo ""
    echo "Usage:"
    echo "./haddjob.sh <job_id>"
    echo "****************************************"
    exit 1
fi

JOB_ID=$1

# Parse the job csv.
while IFS=, read -r id jobname inputfile
do
    if [ "$id" = "$JOB_ID" ]; then

        # Set this to indicate that the job_id entered was actually found.
        # Used purely for the output seen at the end of this file.
        found_job=1
        
        # This is the directory created when running the job_id.
        JOB_DIR="$JOB_DIR_PREFIX$jobname"
        
        # If directory exists
        if [ -d $JOB_DIR ]; then

            # This is the output file created by this script.
            OUTPUT_FILE="output_$jobname.root"

            if [ -f "$OUTPUT_FILE" ]; then
                echo "$OUTPUT_FILE exists already!"
                echo "Make sure the job isn't already done."
            else
                # Move into the directory, hadd the outputs, then move it 
                # one directory above.
                echo "cd $JOB_DIR"
                echo "hadd $OUTPUT_FILE output_${jobname}_*.root"
                echo "mv $OUTPUT_FILE .."
                echo "$OUTPUT_FILE successfully created."
                exit 0
            fi
        else
            echo "Directory '$JOB_DIR' doesn't exist! Make sure to run the job first."
            exit 1
        fi
    fi
done < $JOB_CSV

# If this value is not set, tell the user the job_id was not found.
if [ -z ${found_job+x} ]; then echo "Unable to find job $JOB_ID."; fi
