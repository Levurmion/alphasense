#!/bin/sh

# AF-X6R8D5-F1-model_v4.pdb.gz
# ../json/AF-X6R8D5-F1-predicted_aligned_error_v4.json


cd /bmm/data/alphafold/human/pdb

countp=0
countj=0

for i in AF*.pdb.gz
do
  countp=$(( countp+1 ))
  json="../json/"`echo $i | cut -d- -f3`"-predicted_aligned_error_v4.json"
  if [ ! -f $json ]
  then
    echo Missing JSON for PDB: "$i
    countj=$(( countj+1 ))
  fi

done

echo "Files found: PDB,JSON ="$countp $countj

