
.. _faq:

Hifiasm FAQ
===========


.. contents::
  :local:


How do I get contigs in FASTA?
-------------------------------------
    The FASTA file can be produced from GFA as follows:
    ::

        awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa

Which types of assemblies should I use?
----------------------------------------
    If parental data is available, ``*dip.hap*.p_ctg.gfa`` produced in trio-binning mode should be always preferred. Otherwise if Hi-C data is available, ``*hic.hap*.p_ctg.gfa`` produced in Hi-C mode is the best choice. Both trio-binning mode and Hi-C mode generate fully-phased assemblies. 

    If you only have HiFi reads, hifiasm in default outputs ``*bp.hap*.p_ctg.gfa``. The primary/alternate assemblies can be also produced by using ``--primary``. All these HiFi-only assemblies are not fully-phased. See `blog <https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies>`_ here for more details.

Are inbred/homozygous genomes supported?
--------------------------------------------------------------------------

    Yes, please use the ``-l0`` option to disable purge duplication step.

Are diploid genomes supported?
-------------------------------------
    Yes, most modules of hifiasm are designed for diploid samples, including purge duplication step, partially phased assembly and fully-phased assembly with trio-binning or Hi-C.

Are polyploid genomes supported?
-------------------------------------

    The ``*r_utg.gfa`` and ``*p_utg.gfa`` are lossless so that they also work for polyploid genomes. However, currently the contig-generation modules of hifiasm are designed for diploid samples, which means both the partially phased assembly and the fully-phased assembly does not directly support polyploid genomes. If it is set to >2, the quality of primary assembly for polyploid genomes might be improved. Please use primary assembly for polyploid samples and run multiple rounds of purging steps using third-party tools such as purge_dups. 

Why one Hi-C integrated assembly is larger than another one?
------------------------------------------------------------

    For some samples like human male, the paternal haplotype should be larger than the maternal haplotype. However, if one assembly is much larger than another one, it should be the issues of hifiasm. To fix it, please set smaller value for ``-s`` (default: 0.55). 

    Another possibility is that hifiasm misidentifies coverage threshold for homozygous reads. For instance, hifiasm prints the following information during assembly: 
    ::

        [M::purge_dups] homozygous read coverage threshold: 36

    In this example, hifiasm identifies the coverage threshold for homozygous reads as ``36``. If it is significantly smaller than the homozygous coverage peak, hifiasm will generate two unbalanced assemblies. In this case, please set ``--hom-cov`` to homozygous coverage peak. Please note that tuning ``--hom-cov`` may affect ``*p_utg*gfa`` so that ``*hic*.bin`` should be deleted. Since v0.15.5, hifiasm can detect such changes and renew Hi-C bin files automatically.

    

For Hi-C integrated assembly, why the assembly size of both haplotypes are much larger than the estimated genome size?
------------------------------------------------------------------------------------------------------------------------------
    It is likely that hifiasm misidentifies coverage threshold for homozygous reads. Hifiasm prints the following information for debugging:
    ::

        [M::stat] # heterozygous bases: 645155110; # homozygous bases: 1495396634

    If most bases of a diploid sample are homozygous, the coverage threshold is wrongly determined by hifiasm. For instance, hifiasm prints the following information during assembly: 
    ::

        [M::purge_dups] homozygous read coverage threshold: 36

    In this example, hifiasm identifies the coverage threshold for homozygous reads as ``36``. If it is much smaller than homozygous coverage peak, hifiasm thinks most reads are homozygous and assign them to both assemblies, making both of them much larger than the estimated haploid genome size. In this case, please set ``--hom-cov`` to homozygous coverage peak. Please note that tuning ``--hom-cov`` may affect ``*p_utg*gfa`` so that ``*hic*.bin`` should be deleted. Since v0.15.5, hifiasm can detect such changes and renew Hi-C bin files automatically.


.. _hic-iss:

How can I tweak parameters to improve Hi-C integrated assembly?
---------------------------------------------------------------
    Compared with the HiFi-only assembly or the trio-binning assembly, the Hi-C integrated assembly is a little bit more complex so that you need to take care of the results. See `Why one Hi-C integrated assembly is larger than another one?`_ and `For Hi-C integrated assembly, why the assembly size of both haplotypes are much larger than the estimated genome size?`_ for details on how to fix potential issues.

    There are several other options that may affect the Hi-C integrated assembly. Increasing the values of ``--n-weight``, ``--n-perturb`` and ``--f-perturb`` may improve phasing results but takes longer time. However, tuning ``--l-msjoin`` is tricky. All these options do not affect ``*p_utg*gfa`` so that ``*hic*.bin`` can be reused.

.. _p-large:

Why the size of primary assembly or partially phased assembly is much larger than the estimated genome size? 
---------------------------------------------------------------------------------------------------------------
    It could be because the estimated genome size is incorrect. Another possibility is that hifiasm does not perform enough purging. Setting smaller value for ``-s`` (default: 0.55) or turning ``--hom-cov`` should be helpful. See :ref:`loginter` for more details.


.. _p-hamming:

Why the hamming error rate or the swith error rate of trio-binning assembly is very high?
---------------------------------------------------------------------------------------------------------------
    In rare cases, a potential issue is that a few contigs may misjoin two haplotypes. For example, half of a contig come from mother while another half come from father. Such misjoined contigs can be fixed by manually breaking. The coordinates of problematic regions can be found by A-lines in GFA file or ``yak trioeval -e`` (see `issue 37 <https://github.com/chhylp123/hifiasm/issues/37>`_ for more details). However, if there are many misjoined contigs or the switch/hamming error rate reported by ``yak trioeval`` is very high, users should check if the parental data is correct (see `issue 130 <https://github.com/chhylp123/hifiasm/issues/130#issuecomment-862347943>`_ for more details).

    Another possibility is that there are some unitigs in unitig graph misjoining two haplotypes. Such problematic unitigs might be ignored by the graph-binning strategy. Set smaller value for ``--t-occ`` forcedly remove unitig including unexpected haplotype-specific reads.

Why does hifiasm stuck or crash? 
-------------------------------------
    In most cases, it is caused by the low quality HiFi reads. A good HiFi dataset should have a k-mer plot like `issue10 <https://github.com/chhylp123/hifiasm/issues/10#issuecomment-616213684>`_ or `issue49 <https://github.com/chhylp123/hifiasm/issues/49#issue-729106823>`_. In contrast, low quality HiFi data often lead to weird k-mer plot like `issue93 <https://github.com/chhylp123/hifiasm/issues/93#issue-852259042>`_. Such weird k-mer plots usually indicate insufficient coverage or presence of contaminants. See :ref:`loginter` for more details. If the HiFi data look fine, please raise an issue at the `issue page <https://github.com/chhylp123/hifiasm/issues>`_. 

What's the usage of different bin files in hifiasm?
----------------------------------------------------
    ``*ec.bin``, ``*ovlp.reverse.bin`` and ``*ovlp.source.bin`` save the results of error correction step. ``*hic*bin`` saves the results of Hi-C alignment. Please note that ``*hic*.bin`` should be deleted when tuning any parameters affecting ``*p_utg*gfa``. There are several parameters which does not change ``*p_utg*gfa``, including ``-s``, ``--seed``, ``--n-weight``, ``--n-perturb``, ``--f-perturb`` and ``--l-msjoin``. Since v0.15.5, hifiasm can detect such changes and renew Hi-C bin files automatically. 

Can I generate HiFi-only assembly first, and then add Hi-C or trio data later?
----------------------------------------------------------------------------------------
    Yes, the HiFi-only assembly, Hi-C phased assembly and trio-binning assembly share the same ``*ec.bin``, ``*ovlp.reverse.bin`` and ``*ovlp.source.bin``.

What is the minimum read coverage required for hifiasm?
-------------------------------------------------------
    Usually >=13x HiFi reads per haplotype. Higher coverage might be able to improve the contiguity of assembly.

Why the primary assembly is more contiguous than the fully-phased assemblies and the partially phased assemblies (i.e. ``*.hap*.p_ctg.gfa``)?
----------------------------------------------------------------------------------------------------------------------------------------------------

    For diploid samples, primary assembly usually has greater N50 but at the expense of highly fragmented alternate assembly. From the method view, the primary assembly has an extra joining step, which joins two haplotypes to make primary assembly more contiguous.

    When producing fully-phased assemblies and partially phased assemblies, hifiasm is designed to keep both haplotypes contiguous. It is important for many downstream applications like SV calling.

My assembly is fragmented or not contiguous enough, how do I improve it?
--------------------------------------------------------------------------

    Raising ``-D`` or ``-N`` may improve the resolution of repetitive regions but takes longer time. These two options affect all types of assemblies and usually do not have a negative impact on the assembly quality. In contrast, ``--purge-max`` only affects primary assembly. Setting larger value for ``--purge-max`` makes primary assembly more contiguous but may collapse repeats or segmental duplications.

    If the assembly is too fragmented, users should check if HiFi data is good enough. See `Why does hifiasm stuck or crash?`_ for details. 

How do I avoid misassemblies?
--------------------------------------------------------------------------
    Set smaller value for ``--purge-max``, ``-s`` and ``-O``, or use the ``-u`` option.