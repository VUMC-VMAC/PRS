args<-commandArgs(TRUE)
clumped <- args[1]
sumstats <- args[2]
pvalues <- args[3]

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

##print out the number of variants present at each p value threshold supplied
#split the pvalues from being comma separated
pvalues <- as.numeric(unlist(strsplit(pvalues, ",")))
print(paste("Out of", nrow(summary_snp_overlap), "variants resulting from clumping, there are the following numbers of variants at each p-value threshold:"))
sapply(1:length(pvalues), function(x) paste(sum(as.numeric(summary_snp_overlap$P) <= pvalues[x]), "variants with p values less than or equal to", pvalues[x])) 

#export plink input file
if(grepl("BETA", names(summary_snp_overlap))){
   summary_snp_overlap <- summary_snp_overlap[,c("SNP","A1","BETA", "P")]
} else {
   summary_snp_overlap <- summary_snp_overlap[,c("SNP","A1","OR", "P")]
}
write.table(summary_snp_overlap,paste0(clumped_stem, "_score_input.txt"),quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
