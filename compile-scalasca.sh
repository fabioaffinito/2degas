#!/bin/bash
module purge
module  load autoload scalasca/2.2
cp make.defs.scalasca make.defs
make 
