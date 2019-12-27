#!/bin/bash

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
if [[ $REMOTE == "/"* ]]; then
    REMOTE="$LXPLUS_WORKSPACE$REMOTE"
elif [[ $REMOTE == "." ]]; then
    REMOTE="$LXPLUS_WORKSPACE/"
else
    REMOTE="$LXPLUS_WORKSPACE/$REMOTE"
fi


# The parameters given to the scp command
PARAMS="$LOCAL $LXPLUS_REMOTE:$REMOTE"

echo "****************************************"
echo "(in my best mario voice) HERE WE GO!"

if [ -d "$SOURCE" ]; then
  echo "Copying directory $SOURCE to $DEST"
  scp -r $PARAMS
else
  echo "Copying file $SOURCE to $DEST"
  scp $PARAMS
fi

echo "That's all, folks!"
echo "****************************************"
