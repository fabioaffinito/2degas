#!/bin/bash
module purge
module load autoload intelmpi
make -f makefile.galileo
