#!/bin/bash
module purge
module load intel
module load mkl
cp makedefs/make.defs.omp make.defs
make 
