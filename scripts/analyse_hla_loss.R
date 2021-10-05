# run this script in lohhla env

# initial implementation analysed lohhla results from ~105 samples. Mutation analysis by LOHHLA required >=15 coverage per allele (reduced from default of >=30).
# ~10%? (see below) fail to have even one HLA gene with >=15 coverage per allele.
# Sites that are both mismatching between alleles have >=15 coverage per allele are reported in these files and analysed here.
# setting a low coverage threshold increases the number of samples that can be analysed, but results in reduced certainty of LOH, and therefore fewer positive identifications.

# Methods (from paper):
# A copy number < 0.5, is classified as subject to loss, and thereby indicative of LOH. To avoid over-calling LOH, we also calculate a p value relating to allelic imbalance for each HLA gene. This p value corresponds to the pairwise difference in logR values at mismatch sites between the two HLA homologs, adjusted to ensure each sequencing read is only counted once. Allelic imbalance is determined if p < 0.01 using the paired Student’s t-Test between the two distributions.

# Methods to be used here:
#
# Analysis:
#Tumour purity must be > 30%???
#Alleles must have 5+ mismatch sites.
#Characterise allele loss by Pval_unique, which presumably counts each physical read only once.
#
#To do:
#Stratify data by A,B,C; disease + A,B,C; 
#Then:
#plot loss allele??



# load libraries
library(data.table)
library(tidyverse)

# read tumour type
tumour_type <- fread(cmd = "cat /SAN/colcc/BCI_ICGC/Clinical_Data/specimen.tsv | cut -f 4,22 | egrep '(tumor|submitted_specimen_id)'", sep = "\t", quote = "\"", header = T, stringsAsFactors = F, data.table = F, nThread = 1)
tumour_type <- tumour_type %>% mutate(tumour_type = unlist(lapply(strsplit(tumour_grade, " "), '[', 1)))

# read tumour purity/ploidy estimates
pur_pl <- fread(cmd = "cat /SAN/colcc/alex_work/samples_analysis/*/sequenza/purity_ploidy.txt | grep -v tumorPurity", sep = "\t", quote = "\"", header = F, stringsAsFactors = F, integer64 = "numeric", blank.lines.skip = T, data.table = F, nThread = 1)
colnames(pur_pl) <- c("region", "normalPloidy", "tumorPurity", "tumorPloidy") # assume normalPloidy == 2

# read data
hla_loss <- fread(cmd = "cat /SAN/colcc/alex_work/samples_analysis/*/lohhla/results/*.15.DNA.HLAlossPrediction_CI.xls | grep BAF | sort | uniq && cat /SAN/colcc/alex_work/samples_analysis/*/lohhla/results/*.15.DNA.HLAlossPrediction_CI.xls | grep -v BAF | sort | uniq", sep = "\t", quote = "\"", header = T, stringsAsFactors = F, integer64 = "numeric", blank.lines.skip = T, data.table = F, nThread = 1) # numobs < cn
# note, important to skip blank lines to recognise header
hla_loss <- pur_pl %>%
    mutate(submitted_specimen_id = unlist(lapply(strsplit(region, "_"), function(x) paste(x[1], x[2], sep = "_")))) %>%
    mutate(submitted_specimen_id = gsub("X$", "", submitted_specimen_id)) %>%
    left_join(tumour_type) %>%
    left_join(hla_loss)
system('ps aux | egrep "(MEM|lvande)"') # using ~1.3GB RAM

# a bit of analysis
length(unique(hla_loss$region)) # 180 regions being analysed.
# note, from tracking while running pipeline, 12 samples unable to analyse:
#  4184437X - low map quality
#  4101626, 4116268, 4131213X, 4170686, 4171810 control homozygous at all three loci
#  4108992, 4122063, 4170577, 4175941X, 4177175X, 4197155 T-test error, data are essentially constant (artefact of low coverage)

# how many mismatch sites in each allele?
table(hla_loss$numMisMatchSitesCov) # 0(105) 1(41) 2(30) 3(22) 4(11) 5(10) 6(7) 7(6)... many data points have <5 sites differing between HLA alleles 

# limiting to HLA alleles with 1+ mismatch sites
length(unique(hla_loss$region[which(hla_loss$numMisMatchSitesCov > 0)])) # 100 samples have at least one HLA gene with 1+ sites differing between HLA alleles
length(which(hla_loss$numMisMatchSitesCov > 0)) # in 100 samples, 172 genes have 1+ mismatches between alleles

# limiting to HLA alleles with 5+ mismatch sites
length(unique(hla_loss$region[which(hla_loss$numMisMatchSitesCov >= 5)])) # 54 samples have at least one HLA gene with 5+ mismatches between its alleles
length(which(hla_loss$numMisMatchSitesCov >= 5)) # in 54 samples, 68 genes have 5+ mismatches between alleles

# how many samples experience HLA-LOH? Use 68 alleles in 54 samples and p-value based on support by unique reads supporting an allele
length(unique(hla_loss$region[which(hla_loss$numMisMatchSitesCov >= 5 & hla_loss$PVal_unique < 0.01)])) # 8 samples show some evidence of LOH. note, we haven't thresholded on tumour purity yet.
length(which(hla_loss$numMisMatchSitesCov >= 5 & hla_loss$PVal_unique < 0.01)) # in the 8 samples showing some evidence of LOH, 9 alleles show evidence of LOH. i.e. most have only one allele with sig LOH

# show records with putative allele loss
hla_loss[which(hla_loss$numMisMatchSitesCov >= 5 & hla_loss$PVal_unique < 0.01), c("region", "tumour_type", "tumorPurity", "tumorPloidy", "tumour_grade", "KeptAllele", "LossAllele", "numMisMatchSitesCov", "propSupportiveSites", "PVal_unique")]
#                         region tumour_type tumorPurity tumorPloidy   tumour_grade        KeptAllele        LossAllele numMisMatchSitesCov propSupportiveSites  PVal_unique
#      tumor_4105746_merged_mdup       DLBCL        0.70        2.20 DLBCL grade cb hla_a_02_01_01_01 hla_a_01_01_01_01                  15           100.00000 1.098907e-04
# tumor_4107137_merged_bam_rmdup       DLBCL        0.38        3.00 DLBCL grade cb    hla_b_40_01_02    hla_b_15_27_02                  18           100.00000 4.568868e-03
# tumor_4134005_merged_bam_rmdup          FL        0.79        2.10     FL grade I hla_c_06_02_01_01 hla_c_07_02_01_03                   8           100.00000 8.711317e-03
#     tumor_4139483X_merged_mdup          FL        0.12        2.50     FL grade I hla_c_07_01_01_01 hla_c_06_02_01_01                   5           100.00000 3.919563e-03
#      tumor_4141476_merged_mdup     unknown        0.31        2.10        unknown    hla_a_11_01_01 hla_a_01_01_01_01                  12           100.00000 8.185632e-03
# tumor_4145528_merged_bam_rmdup     unknown        0.47        2.65        unknown    hla_b_07_02_01    hla_b_08_01_01                   5           100.00000 6.942905e-03
# tumor_4178655_merged_bam_rmdup          FL        0.72        2.30     FL grade I    hla_a_30_02_01 hla_a_29_02_01_01                  37           100.00000 4.984888e-08
# tumor_4184094_merged_bam_rmdup       DLBCL        0.80        2.15 DLBCL grade cb hla_a_02_01_01_01 hla_a_24_02_01_01                  35            82.85714 1.215034e-03
# tumor_4184094_merged_bam_rmdup       DLBCL        0.80        2.15 DLBCL grade cb    hla_c_16_02_01 hla_c_05_01_01_02                   8           100.00000 5.966275e-03

# create something to plot
plot_table <- hla_loss %>% filter(numMisMatchSitesCov >= 5) %>%
    mutate(kept_allele_cn = ifelse(KeptAllele == HLA_A_type1, HLA_type1copyNum_withoutBAF, HLA_type2copyNum_withoutBAF)) %>%
    mutate(loss_allele_cn = ifelse(KeptAllele == HLA_A_type1, HLA_type2copyNum_withoutBAF, HLA_type1copyNum_withoutBAF)) %>%
    mutate(hla_gene = substr(HLA_A_type1, 1, 5)) %>%
    select(region, tumour_type, hla_gene, kept_allele_cn, loss_allele_cn, PVal_unique) %>%
    pivot_longer(cols = all_of(c(4,5)), names_to = "which_allele", values_to = "cn") %>%
    group_by(region, hla_gene) %>% mutate(min_cn = min(cn)) %>%
    mutate(sig = ifelse(cn == min_cn & PVal_unique < 0.01, "sig", "ns")) %>% ungroup()
    #%>% mutate(abs_cn = abs(cn)) %>% mutate(max_abs_cn = max(abs_cn)) %>%
    #mutate(sig = ifelse(abs_cn == max_abs_cn & PVal_unique < 0.01, "YES", "NO")) %>% ungroup()

# do summary plot, then extract plots from result pdfs to make overall summary.
plot_table <- plot_table[which(!is.na(plot_table$sig)),]
p <- ggplot(plot_table, aes(x=which_allele, y=cn, col=hla_gene, alpha=sig)) +
    geom_point(width = 0.0, size = 0.25, position = position_jitterdodge(seed = 100, dodge.width = 0.7, jitter.width = 0.2)) + ## size = mm; width = portion of horizontal space to next class
    scale_color_manual(values=c("hla_a" = "#0050a0", "hla_b" = "#800080", "hla_c" = "red4")) +
    scale_alpha_manual(values=c("ns" = 0.2, "sig" = 1)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", colour = "black", size = 0.35, width = 0.2) +
    #stat_summary(aes_string(group = params$colour_var), fun.data = mean_se, geom = "errorbar", size = 0.35, width = 0.2, position = position_dodge(width = 0.75)) + # size = linewidth; width = errorbar width, where each class is at an integer separation; width=1 means adjacent error bars just touch
    stat_summary(fun.y = "mean", geom = "point", pch = 16, colour = "black", size = 0.7) + # , colour = "black"
    #stat_summary(aes_string(group = params$colour_var), fun.y = "mean", geom = "point", pch = 16, colour = "black", size = 0.7, position = position_dodge(width = 0.75)) +  # no 'line' exists as a geom ('line' proceeds across the plot) so use a pch == 3 instead; but this is not centred; use pch 16, but this is offcentre
    theme_classic(base_size = 8) +
    ylab("log CN") +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.box.spacing = unit(1, units = "mm"),
          legend.key.size = unit(2, units = "mm"),
          plot.margin=unit(c(0.15,0.25,0.05,0.05), "cm"),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.title.x = element_blank())
ggsave("~/alex_work/loh_pipeline/scripts/downloaded_results/analyse_hla_loss.pdf", plot = p, width = 55, height = 45, units = "mm", dpi = 300)
#save.image("analyse_hla_loss.rda")


############################################################################################################
#
# below, read CN data. problem: no id to match these with hla-loss records
#
############################################################################################################
cn <- fread(cmd = "cat /SAN/colcc/alex_work/samples_analysis/*/lohhla/results/*.DNA.IntegerCPN_CI.xls | grep BAF | sort | uniq && cat /SAN/colcc/alex_work/samples_analysis/*/lohhla/results/*.DNA.IntegerCPN_CI.xls | grep -v BAF | sort | uniq", sep = "\t", quote = "\"", header = T, stringsAsFactors = F, integer64 = "numeric", blank.lines.skip = T, data.table = F, nThread = 1) # more recs per patient, but no indication of patient id

# can we match log ratios?
length(which(cn$logR_type1 %in% hla_loss$HLAtype1Log2MedianCoverage)) # 11
length(which(cn$logR_type2 %in% hla_loss$HLAtype2Log2MedianCoverage)) # 14
length(which(cn$logR_type1 %in% hla_loss$HLAtype1Log2MedianCoverageAtSites)) # 48
length(which(cn$logR_type2 %in% hla_loss$HLAtype2Log2MedianCoverageAtSites)) # 42
length(which(cn$binlogRtype1 %in% hla_loss$HLAtype1Log2MedianCoverageAtSites)) # 15
length(which(cn$binlogRtype2 %in% hla_loss$HLAtype2Log2MedianCoverageAtSites)) # 17

options(width = 340)
h1 <- hist(hla_loss$HLAtype1Log2MedianCoverageAtSites, breaks = seq(-8, 8, 0.1), plot = F); round(h1$density[50:110], 2)
h2 <- hist(hla_loss$HLAtype2Log2MedianCoverageAtSites, breaks = seq(-8, 8, 0.1), plot = F); round(h2$density[50:110], 2)
dev.new(height = 10, width = 16); matplot(1:length(h1$density), cbind(h1$density, h2$density), type = "l", lty = "solid", lwd = 2, xlim = c(30,100))

data <- hla_loss %>% full_join(cn, by = c("HLAtype1Log2MedianCoverageAtSites" = "logR_type1", "HLAtype2Log2MedianCoverageAtSites" = "logR_type2")) # somewhat ugly; poor matching, if any
data[which(!is.na(data$HLAtype1Log2MedianCoverageAtSites) & !is.na(missMatchseq1)),] # not many good ones

# try looking at individual samples - checking a couple of samples, we find that there is poor or no match between CN analysis and HLA analysis. Also, no indication of allele in the CN analysis.
