#!/bin/bash

if (( $# != 2 )); then
  echo "****************************************"
  echo "Description: SCP a file or directory from local to lxplus7."
  echo ""
  echo "Usage:"
  echo "./brancp.sh <source> <destination>"
  echo ""
  echo "Note: 1. source can be a directory or file."
  echo "      2. destination is within the workspace."
  echo "****************************************"
  exit 1
fi

# this is the username used to connect to lxplus.
NAME="brwillia"

# this is the server incase that changes.
SERVER="lxplus7.cern.ch"

# this be the workspace directory
WORKSPACEDIR="~/nobackup/totem_v2/CMSSW_7_6_3/src/cms-totem-ntuples/Workspace"

# First parameter is the source file (locally)
SOURCE=$1

# Second paramter is the destination file within the above workspace directory.
DEST=$2

# Append the workspace directory to the user-specified directory
# and take into account forward slashes and periods.
if [[ $DEST == "/"* ]]; then
  DEST="$WORKSPACEDIR$DEST"
elif [[ $DEST == "." ]]; then
  DEST="$WORKSPACEDIR/"
else
  DEST="$WORKSPACEDIR/$DEST"
fi

# The parameters given to the scp command
PARAMS="$SOURCE $NAME@$SERVER:$DEST"

echo "(in my best mario voice) HERE WE GO!"

if [ -d "$SOURCE" ]; then
  echo "Copying directory $SOURCE to $DEST"
  scp -r $PARAMS
else
  echo "Copying file $SOURCE to $DEST"
  scp $PARAMS
fi

echo "That's all, folks!"
