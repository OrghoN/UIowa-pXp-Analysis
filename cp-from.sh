#!/bin/bash

if (( $# != 2 )); then
  echo "****************************************"
  echo "Description: SCP a file or directory from lxplus7 to local."
  echo ""
  echo "Usage:"
  echo "./cp-from.sh <remote> <local>"
  echo ""
  echo "Note: 1. remote is within the workspace of lxplus7."
  echo "****************************************"
  exit 1
fi

echo "SCPing from Workspace"

# First parameter is the source file(s) from lxplus7 within the above workspace directory.
REMOTE=$1

# Second parameter is the local destination for the files from lxplus7.
LOCAL=$2

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
PARAMS="$LXPLUS_REMOTE:$REMOTE $LOCAL"
echo $PARAMS

echo "****************************************"

echo "(in my best mario voice) HERE WE GO!"

echo "SCPing from $REMOTE to $LOCAL..."

scp $PARAMS

echo "That's all, folks!"

echo "****************************************"
