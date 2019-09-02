#!/bin/bash
# My first script
if [ $# -eq 1 ]
then
echo "input file is: "$1
echo $1" | paste - - | sort -k1,1 -t \" \" | tr \"\t\" \"\n\" > sort_"$1
cat $1 | paste - - | sort -k1,1 -t " " | tr "\t" "\n" > sort_$1
else
echo "./sort_fasta.sh file_name"
fi
