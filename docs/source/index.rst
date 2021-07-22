Hifiasm
=======

.. toctree::
   :hidden:

   pa-assembly
   trio-assembly
   hic-assembly
   interpreting-output
   faq
   parameter-reference




`Hifiasm <https://github.com/chhylp123/hifiasm>`_ is a fast haplotype-resolved de novo assembler for PacBio HiFi reads. It can assemble a human genome in several hours and assemble a ~30Gb California redwood genome in a few days. Hifiasm emits partially phased assemblies of quality competitive with the best assemblers. Given parental short reads or Hi-C data, it produces arguably the best haplotype-resolved assemblies so far.

Publications
============

Hifiasm
  Haoyu Cheng, Gregory T. Concepcion, Xiaowen Feng, Haowen Zhang & Heng Li.
  `Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm <https://doi.org/10.1038/s41592-020-01056-5>`_. Nature Methods. (2021).

Install
=======
The easiest way to get started is to download a `release <https://github.com/chhylp123/hifiasm/releases>`_. Please report any issues on `github issues <https://github.com/chhylp123/hifiasm/issues>`_ page.

In addition, the latest unreleased version can be found from github:

::

  git clone https://github.com/chhylp123/hifiasm
  cd hifiasm && make

Another way is to install hifiasm via `bioconda <https://anaconda.org/bioconda/hifiasm>`_:

::

  conda install -c bioconda hifiasm

Assembly Concepts 
=================
There are different types of assemblies which are commonly used in practice (see 
`details <https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies>`_). 
Hifiasm produces primary/alternate assemblies or partially phased assemblies 
only with HiFi reads. Given Hi-C data or trio-binning data, hifiasm produces 
contiguous fully-phased assemblies, i.e. haplotype-resolved assemblies.

Why Hifiasm?
============
* Hifiasm delivers high-quality assemblies. It tends to generate longer contigs
  and resolve more segmental duplications than other assemblers.

* Given Hi-C reads or short reads from the parents, hifiasm can produce overall the best
  haplotype-resolved assembly so far. It is the assembler of choice by the 
  `Human Pangenome Project <https://humanpangenome.org/>`_ for the first batch of samples.

* Hifiasm can purge duplications between haplotigs without relying on
  third-party tools such as purge\_dups. Hifiasm does not need polishing tools
  like pilon or racon, either. This simplifies the assembly pipeline and saves
  running time.

* Hifiasm is fast. It can assemble a human genome in half a day and assemble a
  ~30Gb redwood genome in three days. No genome is too large for hifiasm.

* Hifiasm is trivial to install and easy to use. It does not required Python,
  R or C++11 compilers, and can be compiled into a single executable. The
  default setting works well with a variety of genomes.

Learn
=====

*  :ref:`HiFi-only Assembly        <pa-assembly>` - Assembling HiFi reads without additional data types
*  :ref:`Trio-binning Assembly     <trio-assembly>` - Producing fully phased assemblies with HiFi and trio-binning data
*  :ref:`Hi-C Integrated Assembly   <hic-assembly>`   - Producing fully phased assemblies with HiFi and Hi-C data
*  :ref:`Hifiasm Output           <interpreting-output>`   - Interpreting results
*  :ref:`Hifiasm FAQ           <faq>`   - Frequently asked questions
*  :ref:`Hifiasm Parameters <parameter-reference>` - Parameter reference of hifiasm
