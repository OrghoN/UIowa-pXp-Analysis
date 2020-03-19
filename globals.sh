#!/bin/bash
#
#
# These are the global variables used within the shell scripts for running jobs
# within lxplus7.
#
#

# Color codes.
RESET="\e[39m"
GRAY="\e[90m"
GREEN="\e[32m"
YELLOW="\e[93m"
RED="\e[91m"

# Colored prefix used to indicate which command is being executed.
EXECUTING_PREFIX="[${GREEN}EXECUTING${RESET}]:"

# Colored prefix used to indicate information.
INFO_PREFIX="[${YELLOW}INFO${RESET}]:"

# Colored prefix used to indicate an error.
ERROR_PREFIX="[${RED}ERROR${RESET}]:"


# Print out an error message using the defined prefix.
print_err() {
  echo -e "$ERROR_PREFIX $*"
}

# Print out an informative message using the defined prefix.
print_info() {
  echo -e "$INFO_PREFIX $*"
}

# Executes a command and prints a message to inform the user.
execute() {
  echo -e "$EXECUTING_PREFIX $*"
  eval $*
}

# Check if the LXPLUS_USER environment variable has been set.This is useful for
# scripts that are ran locally and depend on user-specific environment 
# variables.
check_vars() {
  if [ -z ${LXPLUS_USER+x} ]; then
    print_err "Environment variables not set!"
    print_err "Make sure to run 'source ./set-tom.sh '"
    print_err "or 'source ./set-bran.sh'"
    exit 1
  fi
}
