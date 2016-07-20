#!/bin/bash
module purge
module load autoload intelmpi
module load mkl
cp makedefs/make.defs.omp make.defs
make 
