#!/bin/bash

source globals.sh
check_vars

execute "ssh -Y $LXPLUS_USER@$LXPLUS_DOMAIN"
