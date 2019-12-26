## Getting Started

```sh
# Install hifiasm (requiring g++ and zlib)
git clone https://github.com/chhylp123/hifiasm.git
cd hifiasm && make
# Assembly (corrected reads are at NA12878.asm.fa, assembly graph is at NA12878.asm.fa.gfa)
./hifiasm -w -l -q NA12878.fq.gz -o NA12878.asm.fa -k 40 -t 32 -r 2
```

## Introduction
Hifiasm is an ultrafast haplotype-resolved de novo assembler based on PacBio Hifi reads. Unlike most existing assemblers, hifiasm starts from uncollapsed genome. Thus, it is able to keep the haplotype information as much as possible. The input of hifiasm is the PacBio Hifi reads in fasta/fastq format, and there are two types of output: (1) haplotype-resolved assembly graph in [GFA][gfa] format; (2) haplotype-aware error corrected reads. So far hifiasm is in early development stage, it will output phased chromosome-level high-quality assembly in the near future.

Hifiasm is a standalone and lightweight assembler, which does not need external libraries (except zlib). For large genomes, it can generate high-quality assembly in a few hours. Hifiasm has been tested on the following datasets:

|<sub>Dataset<sub>|<sub>GSize<sub>|<sub>Cov<sub>|<sub>Asm options<sub>|<sub>CPU time<sub>|<sub>Wall time<sub>|<sub>RAM<sub>|<sub>[unitig][unitig]/contig N50<sup>[1]</sup><sub>|
|:---------------|-----:|-----:|:---------------------|-------:|--------:|----:|----------------:|
|<sub>[Human NA12878]<sub>|<sub>3Gb<sub>|<sub>x28<sub>|<sub>-k 40 -t 42 -r 2<sub>|<sub>200h<sub>|    <sub>5h32m<sub>|<sub>?114G<sub>|<sub>93.5Kb/18.6Mb<sub>|
|<sub>[Human HG002]<sub>|<sub>3Gb<sub>|<sub>x43<sub>|<sub>-k 40 -t 42 -r 2<sub>|<sub>?405h<sub>|<sub>?10h12m<sub>|<sub>?122G<sub>|<sub>320kb/29Mb<sub>|
|<sub>[Human CHM13]<sub>|<sub>3Gb<sub>|<sub>x27<sub>|<sub>-k 40 -t 42 -r 2<sub>|<sub>?<sub>|<sub>?<sub>|<sub>?<sub>|<sub>NA<sup>[2]</sup>/39.8Mb<sub>|
|<sub>[Butterfly]<sub>|<sub>358Mb<sub>|<sub>x35<sub>|<sub>-k 40 -t 42 -r 2 -z 20<sub>|<sub>17.1h<sub>|<sub>36m<sub>|<sub>16G<sub>|<sub>7.5Mb/NA<sup>[3]</sup><sub>|

<sub>[1] unitig N50 is the N50 of assembly graph with haplotype information (i.e., bubbles), while the contig N50 is the N50 of haplotype collapsed assembly (i.e., without bubbles).
[2] CHM13 is a homozygous sample, so unitig N50 makes no sense.
[3] Butterfly has high heterozygous rate, so most chromosomes have been fully separated into two haplotypes. In this case, contig N50 makes no sense.<sub>

## Usage
For Hifi reads assembly, a typical command line looks like:

```sh
./hifiasm -w -l -q NA12878.fq.gz -o NA12878.asm.fa -k 40 -t 32 -r 2
```

where `-q` specifies the input reads and `-o` specifies the output files. In this example, the assembly graph can be found at NA12878.asm.fa.gfa, and the corrected reads can be found at NA12878.asm.fa. `-k`, `-t` and `-r` specify the length of k-mer, the number of CPU threads, and the number of correction rounds, respectively. Note that `-w` means hifiasm will save all overlaps to disk, which can avoid the time-consuming all-to-all overlap calculation next time. For hifiasm with `-l`, if the overlap information has been obtained by `-w` in advance, it is able to load all overlaps from disk and then directly do assembly.

Please note that some old hifi reads may consist of short adapters. To improve the assembly quality, adapters should be removed by `-z` as follow:

```sh
./hifiasm -w -l -q butterfly.fq.gz -o butterfly_asm_clean.fa -k 40 -t 42 -r 2 -z 20
```

In this example, hifiasm will remove 20 bases from both ends of each read.



[unitig]: http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Terminology
[gfa]: https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md