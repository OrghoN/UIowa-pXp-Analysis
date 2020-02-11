#!/bin/bash

# This is the app being submitted.
APPNAME="tomanaRP4k"

# All of the inputs to "submit-condorRP.csh" are located here along with a 
# corresponding job_id.
JOB_CSV="all-trks.csv"

if (( $# != 1 )); then
    echo "****************************************"
    echo "Run one of the jobs for all-trks."
    echo "The job_ids are numbered 1-14."
    echo ""
    echo "Usage:"
    echo "./runjob.sh <job_id>"
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

        echo "sh ./submit-condorRP.csh $APPNAME $jobname $inputfile"
        exit 0
    else
        echo "$id,$jobname,$inputfile"
    fi
done < $JOB_CSV


# If this value is not set, tell the user the job_id was not found.
if [ -z ${found_job+x} ]; then echo "Unable to find job $JOB_ID."; fi
