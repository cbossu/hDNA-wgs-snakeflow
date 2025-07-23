# hDNA-wgs-snakeflow
Howdy folks! Welcome to my snakemake pipeline for processing historical WGS data. The skeleton for this workflow is from [Eriq's mega-non-model-wgs-snakeflow](https://github.com/eriqande/mega-non-model-wgs-snakeflow), but adapted to follow [Sheela's for hDNA processing protocol](https://www.nature.com/articles/s41558-023-01696-3).

## Where I stray from Sheela's protocol
Sheela uses the perl-based prinse-lite program to filter out low-complexity reads prior to mapping. I use the faster C++ implementation of this program (prinseq++), b/c it is faster, but also b/c it is available on conda.

## Finicky bits
SeqPrep2 is not available on conda, so you will have to download the binary for this program and update the [path]() with your own.

I haven't figure out the best way to deal with this problem yet, but depending on the amount of DNA damage in your samples, you may or may not rescale base quality scores. This workflow assumes all samples will need rescaling. If you haved samples with DNA damage levels too low to warrant rescaling, either use the rescaled bams anyways (they shouldn't change much i don't think), or mv the rmdup bams into the rescaled bams folder, touch, and proceed.

You may not be able to call genotypes after mapping with bwa aln + samse. You will need to try to figure out how to move past this.

## Prep configs
Follow this markdown. See the mnm snakeflow for more detail.

## Quick start
```
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake-8.20.4 snakemake=8.20.4
```
