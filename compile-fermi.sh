#!/bin/bash
module purge
module load bgq-xl
module load lapack
module load essl
module load blas
cp make.defs.fermi make.defs
make
