args<-commandArgs(TRUE)
clumped <- args[1]
sumstats <- args[2]

##Script to generate input file for plink --score

library(data.table)

#get the stem
clumped_stem <- gsub("\\..*","",clumped)

#clumped variants
clumped <- fread(clumped, select = "SNP")

#import summary stats
sumstats <- fread(sumstats)

#restrict to list of SNPs from LD clump
summary_snp_overlap <- sumstats[sumstats$SNP %in% clumped$SNP,]

#export plink input file
summary_snp_overlap <- summary_snp_overlap[,c("SNP","A1","BETA")]
write.table(summary_snp_overlap,paste0(clumped_stem, "_score_input.txt"),quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
