
.. _parameter-reference:

Hifiasm Parameter Reference
============================

Synopsis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assembly only with HiFi reads:
::

  hifiasm -o [prefix] -t [nThreads] [options] input1.fq [input2.fq [...]]

Trio binning assembly with yak dumps:
::

  yak count -o paternal.yak -b37 [-t nThreads] [-k kmerLen] paternal.fq.gz
  yak count -o maternal.yak -b37 [-t nThreads] [-k kmerLen] maternal.fq.gz
  hifiasm [-o prefix] [-t nThreads] [options] -1 paternal.yak -2 maternal.yak child.hifi.fq.gz

Hi-C integrated assembly:
::

  hifiasm -o [prefix] -t [nThreads] --h1 [hic_r1.fq.gz,...] --h2 [hic_r2.fq.gz,...] [options] HiFi.read.fq.gz

To get detailed description of options, run:
::

  hifiasm -h

or:
::

  man ./hifiasm.1


General options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _oopt:

**\-o <FILE=hifiasm.asm>**
  Prefix of output files. See :ref:`outfile` and :ref:`outformat` for more details.

.. _topt:

**\-t <INT=1>**
  Number of CPU threads used by hifiasm.

.. _hopt:

**\-h** 
  Show help information.

.. _versionopt:

**\-\-version** 
  Show version number.


Error correction options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _kopt:

**\-k <INT=51>**
  K-mer length. This option must be less than 64.

.. _wopt:

**\-w <INT=51>**    
  Minimizer window size.

.. _fopt:

**\-f <INT=37>**    
  Number  of  bits for bloom filter; 0 to disable. This bloom filter is used to filter out singleton k-mers when counting all k-mers. It takes 2\ :sup:`(INT-3)` bytes of memory. A proper setting saves memory. ``-f37`` is recommended for human assembly. For small genomes, use ``-f0`` to disable the initial bloom filter which takes 16GB memory at the beginning. For genomes much larger than human, applying ``-f38`` or even ``-f39`` is preferred to save memory on k-mer counting.

.. _Dopt:

**\-D <FLOAT=5.0>**    
  Drop k-mers occurring ``>FLOAT*coverage`` times. Hifiasm discards these high-frequency k-mers during error correction to reduce running time. The ``coverage`` is determined automatically by hifiasm based on k-mer plot, representing homozygous read coverage. Raising this option may improve the resolution of repetitive regions but takes longer time.

.. _NEopt:

**\-N <INT=100>**
  Consider up to ``max(-D*coverage,-N)`` overlaps for each oriented read. The ``coverage`` is determined automatically by hifiasm based on k-mer plot, representing homozygous read coverage. Raising this option may improve the resolution of repetitive regions but takes longer time.

.. _ropt:

**\-r <INT=3>**
  Rounds of haplotype-aware error correction. This option affects all outputs of hifiasm. Odd rounds of correction are preferred in practice.


.. _zopt:

**\-z <INT=0>**
  Length  of  adapters that should be removed. This option remove ``INT`` bases from both ends of each read.  Some old HiFi reads may consist of short adapters (e.g. 20bp adapter at one end). For such data, trimming short adapters would significantly improve the assembly quality.

.. _max-kocc-opt:

**\-\-max-kocc <INT=2000>**
  Employ k-mers occurring < ``INT`` times to rescue repetitive overlaps. This option may improve the resolution of repeats.


.. _hg-size-opt:

**\-\-hg-size <INT(k/m/g)>**
  Estimated haploid genome size used for inferring read coverage. This option is used to get accurate homozygous read coverage during error correction. Common suffices are required, for example, 100m or 3g.


.. _min-hist-cnt-opt:

**\-\-min-hist-cnt <INT=5>**
  When analyzing the k-mer spectrum, ignore counts below ``INT``. For very low coverage of HiFi data, set smaller value for this option. See `issue 45 <https://github.com/chhylp123/hifiasm/issues/49>`_ for example.



Assembly options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _aopt:

**\-a <INT=4>**
  Rounds of assembly graph cleaning. This option is used with ``-x`` and ``-y``. Note that unlike -r, this option does not affect error corrected reads and all-to-all overlaps.


.. _mopt:

**\-m <INT=10000000>**
  Maximal probing distance for bubble popping when generating primary/alternate contig graphs. Bubbles longer than ``INT`` bases will not be popped.

.. _popt:

**\-p <INT=0>**
  Maximal probing distance for bubble popping when generating haplotype-resolved processed unitig graph without small bubbles. Bubbles longer than ``INT`` bases will not be popped. Small bubbles might be caused by somatic mutations or noise in data. Please note that hifiasm automatically pops small bubbles based on coverage, which can be tweaked by ``--hom-cov``.

.. _nopt:

**\-n <INT=3>**
  A unitig is considered small if it is composed of less than ``INT`` reads. Hifiasm may try to remove small unitigs at various steps.

.. _xyopt:

**\-x <FLOAT1=0.8>, \-y <FLOAT2=0.2>**
  Max and min overlap drop ratio. This option is used with ``-a``. Given a node N in the assembly graph, let max(N) be the length of the longest overlap of N. Hifiasm iteratively drops overlaps of N if their length/max(N) is below a threshold controlled by ``-x`` and ``-y``. Hifiasm applies ``-a`` rounds of short overlap removal with an increasing threshold between ``FLOAT1`` and ``FLOAT2``.

.. _iopt:

**\-i**
  Ignore all bin files so that hifiasm will start again from scratch.

.. _uopt:

**\-u**
  Disable post-join step for contigs which may improve N50. The post-join step of hifiasm improves contig N50 but may introduce misassemblies. 


.. _hom-cov-opt:

**\-\-hom-cov <INT>**
  Homozygous read coverage inferred automatically in default. This option affects different types of outputs, including Hi-C phased assembly and HiFi-only assembly. For more details, see :ref:`hic-iss`, :ref:`p-large` and :ref:`loginter`.

.. _pri-range-opt:

**\-\-pri-range <INT1[,INT2]>**
  Min and max coverage cutoffs of primary contigs. Keep contigs with coverage in this range at p_ctg.gfa. Inferred automatically in default. If ``INT2`` is not specified, it is set to infinity. Set -1 to disable.

.. _lowQ-opt:

**\-\-lowQ <INT=70>**
  Output contig regions with ``>=INT%`` inconsistency to the bed file with suffix lowQ.bed. Set 0 to disable.

.. _b-cov-opt:

**\-\-b-cov <INT=0>**
  Break contigs at potential misassemblies with ``<INT``-fold coverage. Work with ``--m-rate``. Set 0 to disable.

.. _h-cov-opt:

**\-\-h-cov <INT=-1>**
  Break contigs at potential misassemblies with ``>INT``-fold coverage. Work with ``--m-rate``. Set -1 to disable.

.. _m-rate-opt:

**\-\-m-rate <FLOAT=0.75>**
  Break contigs with ``<=FLOAT*coverage`` exact overlaps. Only work when ``--b-cov`` and ``--h-cov`` are specified.

.. _primary-opt:

**\-\-primary**
  Output a primary assembly and an alternate assembly. Enable this option or ``-l0`` outputs a primary assembly and an alternate assembly.


Trio-binning options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _1opt:

**\-1 <FILE>**
  K-mer dump generated by `yak count <https://github.com/lh3/yak>`_ from the paternal/haplotype1 reads.

.. _2opt:

**\-2 <FILE>**
  K-mer dump generated by `yak count <https://github.com/lh3/yak>`_ from the maternal/haplotype2 reads.

.. _3opt:

**\-3 <FILE>**
  List of paternal/haplotype1 read names.

.. _4opt:

**\-4 <FILE>**
  List of maternal/haplotype2 read names.

.. _cdopt:

**\-c <INT1=2>, -d <INT2=5>**
  Lower bound and upper bound of the binned k-mer's frequency. When doing trio binning, a k-mer is said to be differentiating if it occurs >= ``INT2`` times in one sample but occurs < ``INT1`` times in the other sample.


.. _t-occ-opt:

**\-\-t-occ <INT=60>**
  Forcedly remove unitig including ``>INT`` unexpected haplotype-specific reads without considering graph topology. For more details, see :ref:`p-hamming`.


Purge duplication options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _ldopt:

**\-l <INT=3>**
  Level of purge duplication. 0 to disable, 1 to only purge contained haplotigs, 2 to purge all types of haplotigs, 3 to purge all types of haplotigs in the most aggressive way. In default, 3 for non-trio assembly, 0 for trio-binning assembly. For trio-binning assembly, only level 0 and level 1 are allowed.

.. _sdopt:

**\-s <FLOAT=0.55>**
  Similarity threshold for duplicate haplotigs that should be purged. In default, 0.75 for ``-l1/-l2``, 0.55 for ``-l3``. This option affects both HiFi-only assembly and Hi-C phased assembly. For more details, see :ref:`hic-iss` and :ref:`p-large`.

.. _ovlpdopt:

**\-O <INT=1>**
  Min number of overlapped reads for duplicate haplotigs that should be purged.

.. _purgeopt:

**\-\-purge-max <INT>**
  Coverage upper bound of purge duplication, which is inferred automatically in default. If the coverage of a contig is higher than this bound, don't apply purge duplication. Larger value makes assembly more contiguous but may collapse repeats or segmental duplications.

.. _nhapopt:

**\-\-n\-hap <INT=2>**
  Assumption of haplotype number. If it is set to >2, the quality of primary assembly for polyploid genomes might be improved.



Hi-C integration options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _h1opt:

**\-\-h1 <FILEs>**
  File names of input Hi-C R1 ``[r1_1.fq,r1_2.fq,...]``.

.. _h2opt:

**\-\-h2 <FILEs>**
  File names of input Hi-C R2 ``[r2_1.fq,r2_2.fq,...]``.

.. _n-weightopt:

**\-\-n-weight <INT=3>**
  Rounds of reweighting Hi-C links. Raising this option may improve phasing results but takes longer time.

.. _n-perturbopt:

**\-\-n-perturb <INT=10000>**
  Rounds of perturbation. Increasing this option may improve phasing results but takes longer time.

.. _f-perturbopt:

**\-\-f-perturb <FLOAT=0.1>**
  Fraction to flip for perturbation. Increasing this option may improve phasing results but takes longer time.

.. _seedopt:

**\-\-seed <INT=11>**
  RNG seed.


.. _l-msjoin:

**\-\-l-msjoin <INT=500000>**
  Detect misjoined unitigs of ``>=INT`` in size; 0 to disable.
