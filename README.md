Code for de novo assembly of unmapped reads and calling of ancestral sequences, as shown in [Taliun et al., 2019](https://www.biorxiv.org/content/10.1101/563866v1)


## Requirements

* In `$PATH`:
  * `makecalls_bin/`
  * `assembly_bin/`
  * Python 2.7+
  * Perl 5
  * BWA v0.7.12+
  * SAMtools v1.7+
  * BEDTools v2.25.0+
  * Cutadapt v1.16+
  * liftOver
  * GEM pre-release 3
  * AGE (`age_align` binary, more specifically)
  * ABySS v2.0.2

* UCSC references indexed for BWA-MEM:
  * hg38
  * gorGor5
  * panPan2
  * panTro6
  * ponAbe3

* Reference indexed for GEM-mapper:
  * hg38 with only chromosomes 1-22,X,Y and hard-masked PAR on chromosome Y.

* UCSC chain files:
  * gorGor5ToHg38.over.chain
  * panPan2ToHg38.over.chain
  * panTro6ToHg38.over.chain
  * ponAbe3ToHg38.over.chain

* phiX index for GEM-mapper (provided in `assembly_res/`)