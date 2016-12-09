# Populus Data

### `*.final.vcf`

Downloaded here: http://datadryad.org/resource/doi:10.5061/dryad.0817m

For local ancestry analysis SNPs were called in three chromosomes (6, 12, 15) for 50 reference individuals (25 pure Populus balsamifera and 25 pure P. trichocarpa) and 68 admixed individuals (see Appendix S1, Supporting Information). We sequenced each of the genotypes at an expected coverage ranging from 15X or 30X using the Illumina HiSeq2000. Short reads from the sequencing libraries were independently aligned to the P. trichocarpa version 3 (v3.0) genome using BWA (version 0.6.1-r104) with default parameters. We corrected mate pair metadata and marked duplicate molecules using the FixMateInformation and MarkDuplicates methods in the Picard package (http://picard.sourceforge.net). Reads present in areas surrounding InDels were re-aligned using the IndelRealigner method from GATK (version v1.5-25-gf46f7d0).  Next, we called SNPs and small indels independently using the UnifiedGenotyper method from GATK. SNPs were then filtered to exclude variants within 3bp of any identified variants, having a mapping quality less than 5, and a variant quality less than 30.

### `populus.data.RData`

Generated with this script: https://github.com/BioShock38/Parser/blob/master/R-package/vignettes/VariantAnnotation.Rmd