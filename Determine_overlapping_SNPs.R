args<-commandArgs(TRUE)
sumstats <- args[1] #Path/Name to summary stats file
genotypes <- args[2] #Path/Name of Genetic Data
PRS_tag <- args[3] #outcome

library(data.table)

## compare the summary stats to the input genotypes to get overlapping, non-palindromic variants

# import already cleaned summary stats
summary <- fread(sumstats, select = c("SNP", "BETA", "P"))

# import Cohort SNPs
genetic <- fread(paste0(genotypes,".bim"), select = c(2,5,6))
names(genetic) <- c("SNP","A1","A2")

# determine SNP overlap
summary_genetic_overlap <- merge(genetic, summary, by = "SNP")

# Remove ambiguous SNPs (ie, A/T and G/C SNPs)
ambigSNPs <- summary_genetic_overlap$SNP[(summary_genetic_overlap$A1=="A" & summary_genetic_overlap$A2=="T") | 
	     				 (summary_genetic_overlap$A1=="T" & summary_genetic_overlap$A2=="A") |
					 (summary_genetic_overlap$A1=="G" & summary_genetic_overlap$A2=="C") |
					 (summary_genetic_overlap$A1=="C" & summary_genetic_overlap$A2=="G")]
summary_genetic_overlap <- summary_genetic_overlap[!summary_genetic_overlap$SNP %in% ambigSNPs,]

# export updated summary file
write.table(summary_genetic_overlap,  paste0(PRS_tag,"_summary_stats_updated.txt"),quote=FALSE, row.names=FALSE,  col.names=TRUE, sep=" ")

# export list of summary-genetic overlapping SNPs for plink '--extract'
write.table(summary_genetic_overlap[,c("SNP")], paste(PRS_tag,"_overlapping_SNPs.txt",sep=""),quote=FALSE,  row.names=FALSE, col.names=FALSE, sep=" ")

#print information to the log
print(paste(format(nrow(summary), big.mark = ","), "variants present in", PRS_tag, "summary statistics."))
print(paste(format(nrow(summary_genetic_overlap)+length(ambigSNPs), big.mark = ","), "variants overlap with", genotypes, "."))
print(paste(format(length(ambigSNPs), big.mark = ","), "ambiguous variants removed, leaving", format(nrow(summary_genetic_overlap), big.mark = ","), "for input to PRS calculation."))
