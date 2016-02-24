#!/bin/bash
module purge
module load autoload intelmpi
source $INTEL_HOME/bin/compilervars.sh intel64
make 
