#!/bin/bash 
#
#   tar inputs for 2degas
#
for node in {1..16}
do
   echo "tar for node $node"
   tasks=$(($node*16))
   for (( i=0; i<$tasks; i++))
   do
      #echo "n26.$node.down.tar.$i"
      tar -rf n26.$node.down.tar n26.down.x.$i
      tar -rf n26.$node.up.tar n26.up.x.$i
   done
   gzip n26.$node.down.tar
   gzip n26.$node.up.tar
done
