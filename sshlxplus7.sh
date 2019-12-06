#!/bin/bash

if (( $# != 1)); then
  echo "****************************************"
  echo "Description: Connect to lxplus7 via SSH."
  echo ""
  echo "Usage:"
  echo "./sshlxplus7.sh <username>"
  echo ""
  echo "Note:"
  echo "****************************************"
  exit 1
fi

USERNAME=$1

ssh -Y $USERNAME@lxplus7.cern.ch
