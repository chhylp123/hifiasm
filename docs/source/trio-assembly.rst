
.. _trio-assembly:

Trio-binning Assembly
=====================

When parental short reads are available, hifiasm can also generate a pair of haplotype-resolved assemblies with trio binning. To perform such assembly, you need to count k-mers first with `yak <https://github.com/lh3/yak>`_ and then do assembly::

 yak count -k31 -b37 -t16 -o pat.yak paternal.fq.gz
 yak count -k31 -b37 -t16 -o mat.yak maternal.fq.gz
 hifiasm -o NA12878.asm -t 32 -1 pat.yak -2 mat.yak NA12878.fq.gz

Here ``NA12878.asm.hap1.p_ctg.gfa`` and ``NA12878.asm.hap2.p_ctg.gfa`` give the assemblies for two haplotypes. In the binning mode, hifiasm does not purge haplotig duplicates by default. Because hifiasm reuses saved overlaps, you can generate both primary/alternate assemblies and trio binning assemblies with::

 hifiasm -o NA12878.asm --primary -t 32 NA12878.fq.gz 2> NA12878.asm.pri.log
 hifiasm -o NA12878.asm -t 32 -1 pat.yak -2 mat.yak /dev/null 2> NA12878.asm.trio.log

The second command line will run much faster than the first. The phasing switch error rate and hamming error rate are able to be evaluated quickly by `yak <https://github.com/lh3/yak>`_::

 yak trioeval -t16 pat.yak mat.yak assembly.fa

The W-line and H-line reported by ``yak trioeval`` indicate switch error rate and hamming error rate respectively::

 W       26714   3029448 0.008818
 H       24315   3029885 0.008025

For this example, the switch error rate is 0.8818% and the hamming error rate is 0.8025%.

In rare cases, a potential issue is that a few contigs may misjoin two haplotypes. For example, half of a contig come from mother while another half come from father. Such misjoined contigs can be fixed by 
manually breaking. The coordinates of problematic regions can be found by A-lines in GFA file or ``yak trioeval -e``. However, if there are many misjoined contigs or the switch/hamming error rate reported by ``yak trioeval`` is very high, users should check if the parental data is correct (see `example <https://github.com/chhylp123/hifiasm/issues/130#issuecomment-862347943>`_ for more details).

