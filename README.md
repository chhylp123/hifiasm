## <a name="started"></a>Getting Started

```sh
# Install hifiasm (requiring g++ and zlib)
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make

# Run on test data (use -f0 for small datasets)
wget https://github.com/chhylp123/hifiasm/releases/download/v0.7/chr11-2M.fa.gz
./hifiasm -o test -t4 -f0 chr11-2M.fa.gz 2> test.log
awk '/^S/{print ">"$2;print $3}' test.bp.p_ctg.gfa > test.p_ctg.fa  # get primary contigs in FASTA

# Assemble inbred/homozygous genomes (-l0 disables duplication purging)
hifiasm -o CHM13.asm -t32 -l0 CHM13-HiFi.fa.gz 2> CHM13.asm.log
# Assemble heterozygous genomes with built-in duplication purging
hifiasm -o HG002.asm -t32 HG002-file1.fq.gz HG002-file2.fq.gz

# Hi-C phasing with paired-end short reads in two FASTQ files
hifiasm -o HG002.asm --h1 read1.fq.gz --h2 read2.fq.gz HG002-HiFi.fq.gz

# Trio binning assembly (requiring https://github.com/lh3/yak)
yak count -b37 -t16 -o pat.yak <(cat pat_1.fq.gz pat_2.fq.gz) <(cat pat_1.fq.gz pat_2.fq.gz)
yak count -b37 -t16 -o mat.yak <(cat mat_1.fq.gz mat_2.fq.gz) <(cat mat_1.fq.gz mat_2.fq.gz)
hifiasm -o HG002.asm -t32 -1 pat.yak -2 mat.yak HG002-HiFi.fa.gz
```
See [tutorial][tutorial] for more details. 

## Table of Contents

- [Getting Started](#started)
- [Introduction](#intro)
- [Why Hifiasm?](#why)
- [Usage](#use)
  - [Assembling HiFi reads without additional data types](#hifionly)
  - [Hi-C integration](#hic)
  - [Trio binning](#trio)
  - [Output files](#output)
- [Results](#results)
- [Getting Help](#help)
- [Limitations](#limit)
- [Citing Hifiasm](#cite)

## <a name="intro"></a>Introduction

Hifiasm is a fast haplotype-resolved de novo assembler for PacBio HiFi reads.
It can assemble a human genome in several hours and assemble a ~30Gb California
redwood genome in a few days. Hifiasm emits partially phased assemblies of
quality competitive with the best assemblers. Given parental short reads or
Hi-C data, it produces arguably the best haplotype-resolved assemblies so far.

## <a name="why"></a>Why Hifiasm?

* Hifiasm delivers high-quality assemblies. It tends to generate longer contigs
  and resolve more segmental duplications than other assemblers.

* Given Hi-C reads or short reads from the parents, hifiasm can produce overall the best
  haplotype-resolved assembly so far. It is the assembler of choice by the
  [Human Pangenome Project][hpp] for the first batch of samples.

* Hifiasm can purge duplications between haplotigs without relying on
  third-party tools such as purge\_dups. Hifiasm does not need polishing tools
  like pilon or racon, either. This simplifies the assembly pipeline and saves
  running time.

* Hifiasm is fast. It can assemble a human genome in half a day and assemble a
  ~30Gb redwood genome in three days. No genome is too large for hifiasm.

* Hifiasm is trivial to install and easy to use. It does not required Python,
  R or C++11 compilers, and can be compiled into a single executable. The
  default setting works well with a variety of genomes.

[hpp]: https://humanpangenome.org

## <a name="use"></a>Usage

### <a name="hifionly"></a>Assembling HiFi reads without additional data types

A typical hifiasm command line looks like:
```sh
hifiasm -o NA12878.asm -t 32 NA12878.fq.gz
```
where `NA12878.fq.gz` provides the input reads, `-t` sets the number of CPUs in
use and `-o` specifies the prefix of output files. For this example, the
primary contigs are written to `NA12878.asm.bp.p_ctg.gfa`. 
Since v0.15, hifiasm also produces two sets of
partially phased contigs at `NA12878.asm.bp.hap?.p_ctg.gfa`. This pair of files
can be thought to represent the two haplotypes in a diploid genome, though with
occasional switch errors. The frequency of switches is determined by the
heterozygosity of the input sample.

At the first run, hifiasm saves corrected reads and
overlaps to disk as `NA12878.asm.*.bin`. It reuses the saved results to avoid
the time-consuming all-vs-all overlap calculation next time. You may specify
`-i` to ignore precomputed overlaps and redo overlapping from raw reads.
You can also dump error corrected reads in FASTA and read overlaps in PAF with
```sh
hifiasm -o NA12878.asm -t 32 --write-paf --write-ec /dev/null
```

Hifiasm purges haplotig duplications by default. For inbred or homozygous
genomes, you may disable purging with option `-l0`. Old HiFi reads may contain
short adapter sequences at the ends of reads. You can specify `-z20` to trim
both ends of reads by 20bp. For small genomes, use `-f0` to disable the initial
bloom filter which takes 16GB memory at the beginning. For genomes much larger
than human, applying `-f38` or even `-f39` is preferred to save memory on k-mer
counting.

### <a name="hic"></a>Hi-C integration

Hifiasm can generate a pair of haplotype-resolved assemblies with paired-end
Hi-C reads:
```sh
hifiasm -o NA12878.asm -t32 --h1 read1.fq.gz --h2 read2.fq.gz HiFi-reads.fq.gz
```
In this mode, each contig is supposed to be a haplotig, which by definition
comes from one parental haplotype only. Hifiasm often puts all contigs from the
same parental chromosome in one assembly. It has cleanly separated chrX and
chrY for a human male dataset. Nonetheless, phasing across centromeres is
challenging. Hifiasm is often able to phase entire chromosomes but it may fail 
in rare cases. Also, contigs from different parental chromosomes are randomly mixed as
it is just not possible to phase across chromosomes with Hi-C.

Hifiasm does not perform scaffolding for now. You need to run a standalone
scaffolder such as SALSA or 3D-DNA to scaffold phased haplotigs.

### <a name="trio"></a>Trio binning

When parental short reads are available, hifiasm can also generate a pair of
haplotype-resolved assemblies with trio binning. To perform such assembly, you
need to count k-mers first with [yak][yak] first and then do assembly:
```sh
yak count -k31 -b37 -t16 -o pat.yak paternal.fq.gz
yak count -k31 -b37 -t16 -o mat.yak maternal.fq.gz
hifiasm -o NA12878.asm -t 32 -1 pat.yak -2 mat.yak NA12878.fq.gz
```
Here `NA12878.asm.dip.hap1.p_ctg.gfa` and `NA12878.asm.dip.hap2.p_ctg.gfa` give the two
haplotype assemblies. In the binning mode, hifiasm does not purge haplotig
duplicates by default. Because hifiasm reuses saved overlaps, you can
generate both primary/alternate assemblies and trio binning assemblies with
```sh
hifiasm -o NA12878.asm -t 32 NA12878.fq.gz 2> NA12878.asm.pri.log
hifiasm -o NA12878.asm -t 32 -1 pat.yak -2 mat.yak /dev/null 2> NA12878.asm.trio.log
```
The second command line will run much faster than the first.

### <a name="output"></a>Output files

Hifiasm generates different types of assemblies based on the input data. 
It also writes error corrected reads to the *prefix*.ec.bin binary file and
writes overlaps to *prefix*.ovlp.source.bin and *prefix*.ovlp.reverse.bin.
For more details, please see the complete [documentation][tutorial_output].

## <a name="results"></a>Results

The following table shows the statistics of several hifiasm primary assemblies assembled with v0.12:

|<sub>Dataset<sub>|<sub>Size<sub>|<sub>Cov.<sub>|<sub>Asm options<sub>|<sub>CPU time<sub>|<sub>Wall time<sub>|<sub>RAM<sub>|<sub> N50<sub>|
|:---------------|-----:|-----:|:---------------------|-------:|--------:|----:|----------------:|
|<sub>[Mouse (C57/BL6J)][mouse-data]</sub>|<sub>2.6Gb</sub> |<sub>&times;25</sub>|<sub>-t48 -l0</sub> |<sub>172.9h</sub> |<sub>4.8h</sub> |<sub>76G</sub> |<sub>21.1Mb</sub>|
|<sub>[Maize (B73)][maize-data]</sub>     |<sub>2.2Gb</sub> |<sub>&times;22</sub>|<sub>-t48 -l0</sub> |<sub>203.2h</sub> |<sub>5.1h</sub> |<sub>68G</sub> |<sub>36.7Mb</sub>|
|<sub>[Strawberry][strawberry-data]</sub> |<sub>0.8Gb</sub> |<sub>&times;36</sub>|<sub>-t48 -D10</sub>|<sub>152.7h</sub> |<sub>3.7h</sub> |<sub>91G</sub> |<sub>17.8Mb</sub>|
|<sub>[Frog][frog-data]</sub>             |<sub>9.5Gb</sub> |<sub>&times;29</sub>|<sub>-t48</sub>     |<sub>2834.3h</sub>|<sub>69.0h</sub>|<sub>463G</sub>|<sub>9.3Mb</sub>|
|<sub>[Redwood][redwood-data]</sub>       |<sub>35.6Gb</sub>|<sub>&times;28</sub>|<sub>-t80</sub>     |<sub>3890.3h</sub>|<sub>65.5h</sub>|<sub>699G</sub>|<sub>5.4Mb</sub>|
|<sub>[Human (CHM13)][CHM13-data]</sub>   |<sub>3.1Gb</sub> |<sub>&times;32</sub>|<sub>-t48 -l0</sub> |<sub>310.7h</sub> |<sub>8.2h</sub> |<sub>114G</sub>|<sub>88.9Mb</sub>|
|<sub>[Human (HG00733)][HG00733-data]</sub>|<sub>3.1Gb</sub>|<sub>&times;33</sub>|<sub>-t48</sub>     |<sub>269.1h</sub> |<sub>6.9h</sub> |<sub>135G</sub>|<sub>69.9Mb</sub>|
|<sub>[Human (HG002)][NA24385-data]</sub> |<sub>3.1Gb</sub> |<sub>&times;36</sub>|<sub>-t48</sub>     |<sub>305.4h</sub> |<sub>7.7h</sub> |<sub>137G</sub>|<sub>98.7Mb</sub>|

[mouse-data]:      https://www.ncbi.nlm.nih.gov/sra/?term=SRR11606870
[maize-data]:      https://www.ncbi.nlm.nih.gov/sra/?term=SRR11606869
[strawberry-data]: https://www.ncbi.nlm.nih.gov/sra/?term=SRR11606867
[frog-data]:       https://www.ncbi.nlm.nih.gov/sra?term=(SRR11606868)%20OR%20SRR12048570
[redwood-data]:    https://www.ncbi.nlm.nih.gov/sra/?term=SRP251156
[CHM13-data]:      https://www.ncbi.nlm.nih.gov/sra?term=(((SRR11292120)%20OR%20SRR11292121)%20OR%20SRR11292122)%20OR%20SRR11292123

Hifiasm can assemble a 3.1Gb human genome in several hours or a ~30Gb hexaploid
redwood genome in a few days on a single machine. For trio binning assembly:

|<sub>Dataset<sub>|<sub>Cov.<sub>|<sub>CPU time<sub>|<sub>Elapsed time<sub>|<sub>RAM<sub>|<sub> N50<sub>|
|:---------------|-----:|-------:|--------:|----:|----------------:|
|<sub>[HG00733][HG00733-data], [\[father\]][HG00731-data], [\[mother\]][HG00732-data]</sub>|<sub>&times;33</sub>|<sub>269.1h</sub>|<sub>6.9h</sub>|<sub>135G</sub>|<sub>35.1Mb (paternal), 34.9Mb (maternal)</sub>|
|<sub>[HG002][NA24385-data],   [\[father\]][NA24149-data], [\[mother\]][NA24143-data]</sup>|<sub>&times;36</sub>|<sub>305.4h</sub>|<sub>7.7h</sub>|<sub>137G</sub>|<sub>41.0Mb (paternal), 40.8Mb (maternal)</sub>|

<!--
|<sub>[NA12878][NA12878-data], [\[father\]][NA12891-data], [\[mother\]][NA12892-data]</sub>|<sub>&times;30</sub>|<sub>180.8h</sub>|<sub>4.9h</sub>|<sub>123G</sub>|<sub>27.7Mb (paternal), 27.0Mb (maternal)</sub>|
-->

[HG00733-data]: https://www.ebi.ac.uk/ena/data/view/ERX3831682
[HG00731-data]: https://www.ebi.ac.uk/ena/data/view/ERR3241754
[HG00732-data]: https://www.ebi.ac.uk/ena/data/view/ERR3241755
[NA24385-data]: https://www.ncbi.nlm.nih.gov/sra?term=(((SRR10382244)%20OR%20SRR10382245)%20OR%20SRR10382248)%20OR%20SRR10382249
[NA24149-data]: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/HG003Run01-13262252/
[NA24143-data]: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004Run01-15133132/
[NA12878-data]: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/PacBio_SequelII_CCS_11kb/
[NA12891-data]: https://www.ebi.ac.uk/ena/data/view/ERR194160
[NA12892-data]: https://www.ebi.ac.uk/ena/data/view/ERR194161

Human assemblies above can be acquired [from Zenodo][zenodo-human] and
non-human ones are available [here][zenodo-nonh].

[zenodo-human]: https://zenodo.org/record/4393631
[zenodo-nonh]: https://zenodo.org/record/4393750
[unitig]: http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Terminology
[gfa]: https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md
[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[yak]: https://github.com/lh3/yak
[tutorial]: https://hifiasm.readthedocs.io/en/latest/index.html
[tutorial_output]: https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output


## <a name="help"></a>Getting Help

For detailed description of options, please see [tutorial][tutorial] or `man ./hifiasm.1`. The `-h`
option of hifiasm also provides brief description of options. If you have
further questions, please raise an issue at the [issue
page](https://github.com/chhylp123/hifiasm/issues).

## <a name="limit"></a>Limitations

1. Purging haplotig duplications may introduce misassemblies.

## <a name="cite"></a>Citating Hifiasm

If you use hifiasm in your work, please cite:

> Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., Li H. (2021)
> Haplotype-resolved de novo assembly using phased assembly graphs with
> hifiasm. *Nat Methods*, **18**:170-175.
> https://doi.org/10.1038/s41592-020-01056-5
