
.. _hic-assembly:

Hi-C Integrated Assembly
========================

Hifiasm can generate a pair of haplotype-resolved assemblies with paired-end Hi-C reads::

 hifiasm -o NA12878.asm -t32 --h1 read1.fq.gz --h2 read2.fq.gz HiFi-reads.fq.gz

In this mode, each contig is supposed to be a haplotig, which by definition comes from one parental haplotype only. Hifiasm often puts all contigs from the same parental chromosome in one assembly. It has cleanly separated chrX and chrY for a human male dataset. Nonetheless, phasing across centromeres is challenging. Hifiasm is often able to phase entire chromosomes but it may fail in rare cases. Also, contigs from different parental chromosomes are randomly mixed as it is just not possible to phase across chromosomes with Hi-C. Hifiasm does not perform scaffolding for now. You need to run a standalone scaffolder such as SALSA or 3D-DNA to scaffold phased haplotigs. 


For samples with high heterozygosity rate, a common issue is that one assembly is much larger than another one. To fix this issue, please set smaller value for ``-s`` (default: 0.55). Another possibility is that hifiasm misidentifies coverage threshold for homozygous reads. In this case, please set ``--hom-cov`` to homozygous coverage peak. See :ref:`hic-iss` for more details.

At the first run, hifiasm saves the alignment of Hi-C reads to disk as ``*hic*.bin``. It reuses the saved results to avoid Hi-C alignment next time. Please note that ``*hic*.bin`` should be deleted when tuning any parameters affecting ``*p_utg*gfa``. Since v0.15.5, hifiasm can detect such changes and renew Hi-C bin files automatically. There are several parameters which do not change ``*p_utg*gfa``, including ``-s``, ``--seed``, ``--n-weight``, ``--n-perturb``, ``--f-perturb`` and ``--l-msjoin``.
