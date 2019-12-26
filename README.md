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

|Dataset         |GSize |  Cov |          Asm options |CPU time|Wall time|  RAM|unitig/contig N50<sup>[1]</sup>|
|:---------------|-----:|-----:|:---------------------|-------:|--------:|----:|----------------:|
|[Human NA12878] |   3Gb|x28   |-k 40 -t 42 -r 2      |    200h|    5h32m|?114G|    93.5Kb/18.6Mb|
|[Human HG002]   |   3Gb|x43   |-k 40 -t 42 -r 2      |   ?405h|  ?10h12m|?122G|       320kb/29Mb|
|[Human CHM13]   |   3Gb|x27   |-k 40 -t 42 -r 2      |       ?|        ?|    ?|         NA<sup>[2]</sup>/39.8Mb|
|[Butterfly]     | 358Mb|x35   |-k 40 -t 42 -r 2 -z 20|   17.1h|      36m|  16G|         NA<sup>[3]</sup>/7.5Mb|


<style>
  .markdown-body table td {
    font-size: 1px !important;
  }
</style>

