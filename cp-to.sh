#!/bin/bash

source globals.sh
check_vars

if (( $# != 2 )); then
  echo "****************************************"
  echo "Description: SCP a file or directory from local to lxplus7."
  echo ""
  echo "Usage:"
  echo "./cp-to.sh <local> <remote>"
  echo ""
  echo "Note: 1. local can be a directory or file."
  echo "      2. remote is within the workspace of lxplus7."
  echo "****************************************"
  exit 1
fi

# First parameter is the source file (locally)
LOCAL=$1

# Second paramter is the destination file within the above workspace directory.
REMOTE=$2

# Append the workspace directory to the user-specified directory
# and take into account forward slashes and periods.

if [[ $REMOTE == "/"* ]]; then # True if the REMOTE begins with a '/'.
    REMOTE="$LXPLUS_WORKSPACE$REMOTE"
elif [[ $REMOTE == "." ]]; then # True if the REMOTE is a single '.'
    REMOTE="$LXPLUS_WORKSPACE/"
else
    REMOTE="$LXPLUS_WORKSPACE/$REMOTE"
fi


# The parameters given to the scp command
PARAMS="$LOCAL $LXPLUS_REMOTE:$REMOTE"

if [ -d "$SOURCE" ]; then
  print_info "Copying directory '$LOCAL' to '$REMOTE' ..."
  execute "scp -r $PARAMS"
else
  print_info "Copying file '$LOCAL' to '$REMOTE'..."
  execute "scp $PARAMS"
fi

