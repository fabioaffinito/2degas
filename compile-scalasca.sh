#!/bin/bash
module purge
module load intel parastation/intel-5.1.4-1_1_g064e3f7
module load UNITE
module load scalasca/2.0b3-mpich2-intel
cp make.defs.scalasca make.defs
make 

