#!/bin/bash 
#
# copy restart files to each SDV node
# nodes available from PBS
#
echo "Copying restart files for nvme dir"
hosts=$(uniq ${PBS_NODEFILE})
for host in $hosts
do
   echo $host
   ssh $host mkdir -p /nvme/tmp/$USER
   scp n26.down.x.* n26.up.x.*  $USER@$host:/nvme/tmp/$USER
#   scp n26.$nodes.tar.gz $USER@$host:/nvme/tmp/$USER
done
rm n26.down.x.* n26.up.x.*

