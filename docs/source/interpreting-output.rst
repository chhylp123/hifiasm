
.. _interpreting-output:

Hifiasm Output 
===============

.. _outfile:

Output files
---------------------------------------

In general, hifiasm generates the following assembly graphs in the GFA format:

* ```prefix`.r_utg.gfa``: haplotype-resolved raw unitig graph. This graph keeps all haplotype information.
* ```prefix`.p_utg.gfa``: haplotype-resolved  processed  unitig  graph  without small bubbles. Small bubbles might be caused by somatic mutations or noise in data, which are not the real haplotype information. Hifiasm automatically pops such small bubbles based on coverage. The option ``--hom-cov`` affects the result. See :ref:`homozygous coverage setting <homcov>` for more details. In addition, the option ``-p`` forcedly pops bubbles.
* ```prefix`.p_ctg.gfa``: assembly graph of primary contigs. This graph includes a complete assembly with long stretches of phased blocks.
* ```prefix`.a_ctg.gfa``: assembly graph of alternate contigs. This graph consists of all contigs that are discarded in primary contig graph.
* ```prefix`.*hap*.p_ctg.gfa``: phased contig graph. This graph keeps the phased contigs.


Hifiasm outputs ``*.r_utg.gfa`` and ``*.p_utg.gfa`` in any cases. Specifically, hifiasm outputs the following assembly graphs in trio-binning mode:

* ```prefix`.dip.hap1.p_ctg.gfa``: fully phased paternal/haplotype1 contig graph keeping the phased paternal/haplotype1 assembly.
* ```prefix`.dip.hap2.p_ctg.gfa``: fully phased maternal/haplotype2 contig graph keeping the phased maternal/haplotype2 assembly.

With Hi-C partition options, hifiasm outputs:

* ```prefix`.hic.p_ctg.gfa``: assembly graph of primary contigs.
* ```prefix`.hic.hap1.p_ctg.gfa``: fully phased contig graph of haplotype1 where each contig is fully phased.
* ```prefix`.hic.hap2.p_ctg.gfa``: fully phased contig graph of haplotype2 where each contig is fully phased.
* ```prefix`.hic.a_ctg.gfa`` (optional with ``--primary``): assembly graph of alternate contigs.

Hifiasm generates the following assembly graphs only with HiFi reads in default:

* ```prefix`.bp.p_ctg.gfa``: assembly graph of primary contigs.
* ```prefix`.bp.hap1.p_ctg.gfa``: partially phased contig graph of haplotype1.
* ```prefix`.bp.hap2.p_ctg.gfa``: partially phased contig graph of haplotype2.

If the option ``--primary`` or ``-l0`` is specified, hifiasm outputs:

* ```prefix`.p_ctg.gfa``: assembly graph of primary contigs.
* ```prefix`.a_ctg.gfa``: assembly graph of alternate contigs.

For each graph, hifiasm also outputs a simplified version (``*noseq*gfa``) without sequences for the ease of visualization. The coordinates of low quality regions are written to ``*lowQ.bed`` in BED format.
The concepts of different types of assemblies can be found `here <https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies>`_.

.. _outformat:

Output file formats
---------------------------------------
Hifiasm broadly follows the specification for `GFA 1.0 <https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md>`_. There are several fields that are specifically used by hifiasm. For ``S`` segment line:

* ``rd:i:``: read coverage. It is calculated by the reads coming from the same contig/unitig.

Hifiasm outputs ``A`` lines including the information of reads which are used to construct contig/unitig. Each ``A`` line is plain-text, tab-separated, and the columns appear in the following order:

.. list-table::
   :widths: 10 25 50
   :header-rows: 1

   * - Col
     - Type
     - Description
   * - 1
     - string
     - Should be always ``A``
   * - 2
     - string
     - Contig/unitig name
   * - 3
     - int
     - Contig/unitig start coordinate of subregion constructed by read
   * - 4
     - char
     - Read strand: "+" or "-"
   * - 5
     - string
     - Read name
   * - 6
     - int
     - Read start coordinate of subregion which is used to construct contig/unitig
   * - 7
     - int
     - Read end coordinate of subregion which is used to construct contig/unitig
   * - 8
     - id:i:int
     - Read ID
   * - 9
     - HG:A:char
     - Haplotype status of read. ``HG:A:a``, ``HG:A:p``, ``HG:A:m`` indicate read is non-binnable, father/hap1-specific and mother/hap2-specific, respectively.

.. _loginter:

Hifiasm log interpretation
---------------------------------------
Hifiasm prints several information for quick debugging, including:

.. _homcov:

* k-mer plot: showing how many k-mers appear a certain number of times. For homozygous samples, there should be one peak around read coverage. For heterozygous samples, there should two peaks, where the smaller peak is around the heterozygous read coverage and the larger peak is around the homozygous read coverage. For example, `issue10 <https://github.com/chhylp123/hifiasm/issues/10#issuecomment-616213684>`_ indicates the heterozygous read coverage and the homozygous read coverage are 28 and 57, respectively. `Issue49 <https://github.com/chhylp123/hifiasm/issues/49#issue-729106823>`_ is another good example. Weird k-mer plot like `issue93 <https://github.com/chhylp123/hifiasm/issues/93#issue-852259042>`_ is often caused by insufficient coverage or presence of contaminants. 
* homozygous coverage: coverage threshold for homozygous reads. Hifiasm prints it as: ``[M::purge_dups] homozygous read coverage threshold: X``. If it is not around homozygous coverage, the final assembly might be either too large or too small. To fix this issue, please set ``--hom-cov`` to homozygous coverage.  
* number of het/hom bases: how many bases in unitig graph are heterozygous and homozygous during Hi-C phased assembly. Hifiasm prints it as: ``[M::stat] # heterozygous bases: X; # homozygous bases: Y``. Given a heterozygous sample, if there are much more homozygous bases than heterozygous bases, hifiasm fails to identify correct coverage threshold for homozygous reads. In this case, please set ``--hom-cov`` to homozygous coverage. 
