#!/bin/bash
module purge
module list
module load autoload intelmpi
cp make.defs.galileo make.defs
make 
