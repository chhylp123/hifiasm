## Getting Started

```sh
# Install hifiasm (requiring g++ and zlib)
git clone https://github.com/chhylp123/hifiasm.git
cd hifiasm && make
# Assembly (corrected reads are at NA12878.asm.fa, assembly graph is at NA12878.asm.fa.gfa)
./hifiasm -w -l -q NA12878.fq.gz -o NA12878.asm.fa -k 40 -t 32 -r 2
```

## Introduction
Hifiasm is an ultrafast haplotype-resolved de novo assembler based on PacBio Hifi reads. Unlike most existing assemblers, hifiasm starts from uncollapsed genome. Thus, it is able to keep the haplotype information as much as possible. The input of hifiasm is the PacBio Hifi reads in fasta/fastq format, and there are two types of output: (1) haplotype-resolved assembly graph in GFA format; (2) haplotype-aware error corrected reads. So far hifiasm is in early development stage, it will output phased chromosome-level high-quality assembly in the near future.

Hifiasm is a standalone and lightweight assembler, which does not need external libraries (except zlib). For large genomes, it can generate high-quality assembly in a few hours. Hifiasm has been tested on the following datasets:
