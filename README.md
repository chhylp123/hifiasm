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

|<sub>Dataset<sub>|<sub>GSize<sub>|<sub>Cov<sub>|<sub>Asm options<sub>|<sub>CPU time<sub>|<sub>Wall time<sub>|<sub>RAM<sub>|<sub>unitig/contig N50<sup>[1]</sup><sub>|
|:---------------|-----:|-----:|:---------------------|-------:|--------:|----:|----------------:|
|<sub>[Human NA12878]<sub>|<sub>3Gb<sub>|<sub>x28<sub>|<sub>-k 40 -t 42 -r 2<sub>|<sub>200h<sub>|    <sub>5h32m<sub>|<sub>?114G<sub>|<sub>93.5Kb/18.6Mb<sub>|
|<sub>[Human HG002]<sub>|<sub>3Gb<sub>|<sub>x43<sub>|<sub>-k 40 -t 42 -r 2<sub>|<sub>?405h<sub>|<sub>?10h12m<sub>|<sub>?122G<sub>|<sub>320kb/29Mb<sub>|
|<sub>[Human CHM13]<sub>|<sub>3Gb<sub>|<sub>x27<sub>|<sub>-k 40 -t 42 -r 2<sub>|<sub>?<sub>|<sub>?<sub>|<sub>?<sub>|<sub>NA<sup>[2]</sup>/39.8Mb<sub>|
|<sub>[Butterfly]<sub>|<sub>358Mb<sub>|<sub>x35<sub>|<sub>-k 40 -t 42 -r 2 -z 20<sub>|<sub>17.1h<sub>|<sub>36m<sub>|<sub>16G<sub>|<sub>NA<sup>[3]</sup>/7.5Mb<sub>|


<style>
  .markdown-body table td {
    font-size: 1px !important;
  }
</style>

