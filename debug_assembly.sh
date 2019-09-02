#!/bin/bash
# My first script
if [ $# -eq 1 ]
then
echo "../minimap2/minimap2 -ax asm20 -t 32 ../minimap2/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz "$1" >"$1".sam"
../minimap2/minimap2 -ax asm20 -t 32 ../minimap2/Homo_sapiens.GRCh38.dna.primary_assembly.fa $1 >$1.sam
echo "samtools view -Sb "$1".sam >"$1".bam"
samtools view -Sb $1.sam >$1.bam
echo "samtools sort "$1".bam sort_"$1
samtools sort $1.bam sort_$1
echo "rm sort_"$1".bam.bai"
rm sort_$1.bam.bai
echo "samtools index sort_"$1".bam"
samtools index sort_$1.bam
else
echo "debug_assembly.sh intput.fa"
fi
