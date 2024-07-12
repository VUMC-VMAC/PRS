args<-commandArgs(TRUE)
sumstats <- args[1] #Path/Name to summary stats file
genotypes <- args[2] #Path/Name of Genetic Data
PRS_tag <- args[3] #outcome
logistic <- args[4] # yes or no indicating whether test stat should be BETA or OR

library(data.table)

setDTthreads(4)

## compare the summary stats to the input genotypes to get overlapping, non-palindromic variants

# import already cleaned summary stats
if(logistic == "no"){
     summary <- fread(sumstats, select = c("SNP", "A1", "A2", "BETA", "P"))
} else if (logistic == "yes"){
     summary <- fread(sumstats, select = c("SNP", "A1", "A2", "OR", "P"))
}

# import Cohort SNPs
genetic <- fread(paste0(genotypes,".bim"), select = c(2,5,6))
names(genetic) <- c("SNP", "A1", "A2")

# Remove ambiguous SNPs (ie, A/T and G/C SNPs)
ambigSNPs <- genetic$SNP[(genetic$A1=="A" & genetic$A2=="T") | 
	     		 (genetic$A1=="T" & genetic$A2=="A") |
			 (genetic$A1=="G" & genetic$A2=="C") |
			 (genetic$A1=="C" & genetic$A2=="G")]
genetic <- genetic[!genetic$SNP %in% ambigSNPs,]

# determine overlap by merging on variant id, A1 and A2
summary_genetic_overlap <- merge(summary, genetic, by = c("SNP", "A1", "A2"))

#check for variants which are just on the opposite strand
summary_genetic_overlap_oppallele <- merge(summary, genetic, by.x = c("SNP", "A1", "A2"), by.y = c("SNP", "A2", "A1"))
summary_genetic_overlap <- rbind(summary_genetic_overlap, summary_genetic_overlap_oppallele)
rm(summary_genetic_overlap_oppallele)

#check for variants which have to be flipped
genetic$A1_old <- genetic$A1
genetic$A2_old <- genetic$A2
genetic$A1[genetic$A1_old == "A"] <- "T"
genetic$A1[genetic$A1_old == "T"] <- "A"
genetic$A1[genetic$A1_old == "G"] <- "C"
genetic$A1[genetic$A1_old == "C"] <- "G"
genetic$A2[genetic$A2_old == "A"] <- "T"
genetic$A2[genetic$A2_old == "T"] <- "A"
genetic$A2[genetic$A2_old == "G"] <- "C"
genetic$A2[genetic$A2_old == "C"] <- "G"
genetic$A1_old <- NULL
genetic$A2_old <- NULL
summary_genetic_overlap_flip <- merge(summary[!summary$SNP %in% summary_genetic_overlap$SNP,], genetic, by = c("SNP", "A1", "A2"))
summary_genetic_overlap_flip_oppallele <- merge(summary, genetic, by.x = c("SNP", "A1", "A2"), by.y = c("SNP", "A2", "A1"))
summary_genetic_overlap_flip <- rbind(summary_genetic_overlap_flip, summary_genetic_overlap_flip_oppallele)
rm(summary_genetic_overlap_flip_oppallele)

#write out the list of those variants to flip at the front end
write.table(summary_genetic_overlap_flip$SNP,  paste(PRS_tag,"_SNPs_to_flip.txt",sep=""),quote=FALSE,  row.names=FALSE, col.names=FALSE, sep=" ")

#combine to keep all the variants which overlap and have matching (or flipped) alleles
summary_genetic_overlap <- rbind(summary_genetic_overlap, summary_genetic_overlap_flip)

# export updated summary file
write.table(summary_genetic_overlap,  paste0(PRS_tag,"_summary_stats_updated.txt"),quote=FALSE, row.names=FALSE,  col.names=TRUE, sep=" ")

# export list of summary-genetic overlapping SNPs for plink '--extract'
write.table(summary_genetic_overlap[,c("SNP")], paste(PRS_tag,"_overlapping_SNPs.txt",sep=""),quote=FALSE,  row.names=FALSE, col.names=FALSE, sep=" ")

#print information to the log
print(paste(format(nrow(summary), big.mark = ","), "variants present in", PRS_tag, "summary statistics."))
print(paste(format(nrow(summary_genetic_overlap)+length(ambigSNPs), big.mark = ","), "variants overlap with", genotypes, ", including", format(nrow(summary_genetic_overlap_flip), big.mark = ","), "variants which need alleles flipped."))
print(paste(format(length(ambigSNPs), big.mark = ","), "ambiguous variants removed, leaving", format(nrow(summary_genetic_overlap), big.mark = ","), "for input to PRS calculation."))
