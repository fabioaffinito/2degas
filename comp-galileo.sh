#!/bin/bash
module purge
module load autoload intelmpi
module load mkl
cp makedefs/make.defs.galileo make.defs
make 
