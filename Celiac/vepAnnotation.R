################################################################################
library(labnotebook20162017)

################################################################################
# load data

load("../../Celiac/lfmm.bridge.res.RData")

################################################################################
# run of vep

chm6.vcf <- lfmm.bridge.res$locus.covariable %>%
  transmute(chromosome = chromosome,
            start = physical.pos,
            end = physical.pos,
            allele = paste0(allele1,"/",allele2),
            strand = NA,
            identifier = marker.ID)
chm6.vep <- variant_effect_predictor(chm6.vcf)
save(chm6.vep, file = "chm6.vep.RData")
