
.. _pa-assembly:

HiFi-only Assembly
==================

A typical hifiasm command line looks like::

 hifiasm -o NA12878.asm -t 32 NA12878.fq.gz

where ``NA12878.fq.gz`` provides the input reads, ``-t`` sets the number of CPUs in
use and ``-o`` specifies the prefix of output files. Input sequences should be FASTA 
or FASTQ format, uncompressed or compressed with gzip (.gz). The quality scores of reads 
in FASTQ are ignored by hifiasm. Hifiasm outputs assemblies in `GFA <https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md>`_ format.

At the first run, hifiasm saves corrected reads and overlaps to disk as ``NA12878.asm.*.bin``. It reuses the saved results to avoid the time-consuming all-vs-all overlap calculation next time. You may specify ``-i`` to ignore precomputed overlaps and redo overlapping from raw reads. You can also dump error corrected reads in FASTA and read overlaps in PAF with::

 hifiasm -o NA12878.asm -t 32 --write-paf --write-ec /dev/null

Hifiasm purges haplotig duplications by default. For inbred or homozygous genomes, you may disable purging with option ``-l0``. Old HiFi reads may contain short adapter sequences at the ends of reads. You can specify ``-z20`` to trim both ends of reads by 20bp. For small genomes, use ``-f0`` to disable the initial bloom filter which takes 16GB memory at the beginning. For genomes much larger than human, applying ``-f38`` or even ``-f39`` is preferred to save memory on k-mer counting.


Produce two partially phased assemblies
---------------------------------------


Since v0.15, hifiasm produces two sets of partially phased contigs in default like::

 hifiasm -o NA12878.asm -t 32 NA12878.fq.gz

In this example, the partially phased contigs are written to ``NA12878.asm.bp.hap*.p_ctg.gfa``. 
This pair of files can be thought to represent the two haplotypes in a diploid genome, though with occasional switch errors. The frequency of switches is determined by the heterozygosity of the input sample. Hifiasm also writes the primary contigs to ``NA12878.asm.bp.p_ctg.gfa``. 

For samples with high heterozygosity rate, a common issue is that one set of partially phased contigs is much larger than another set. To fix this issue, please set smaller value for ``-s`` (default: 0.55). Another possibility is that hifiasm misidentifies coverage threshold for homozygous reads. 
In this case, please set ``--hom-cov`` to homozygous coverage. See :ref:`p-large` for more details.


Produce primary/alternate assemblies
------------------------------------

To get primary/alternate assemblies, the option ``--primary`` should be set::

 hifiasm -o NA12878.asm --primary -t 32 NA12878.fq.gz

The primary contigs and the alternate contigs are written to ``NA12878.asm.p_ctg.gfa`` and ``NA12878.asm.a_ctg.gfa``, respectively. For inbred or homozygous genomes, the primary/alternate assemblies can be also produced by ``-l0``. Similarly, turning ``-s`` or ``--hom-cov`` should 
be helpful if the primary assembly is too large. See :ref:`p-large` for more details.

