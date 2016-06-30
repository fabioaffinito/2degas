#!/bin/bash 
#
# copy restart files to each SDV node
# nodes available from PBS
#
echo "Removing files from /nvme directory"
hosts=$(uniq ${PBS_NODEFILE})
for host in $hosts
do
  echo $host
  ssh $host rm /nvme/tmp/$USER/*
done

