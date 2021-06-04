#cleaning kunkle summary stats

library(data.table)

setDTthreads(4)

##read in summary stats, n variants = 9,556,884
sumstats <- fread("/nfs/clarklc/PRS/Summary_Stats/IGAP2_ModelB_METAL_COMMON.InvVar.woADNI_BIOCARD_WASHU.METAANALYSIS.TBL")

##remove multi-allelic variants
#split into chr:bp to be able to detect multi-allelic variants
sumstats$chrbp <- sub(':([^:]*)$', '',sumstats$MarkerName)

#remove variants at the same position, n = 9,535,830
num_multiallelic <- sum(duplicated(sumstats$chrbp))
sumstats <- sumstats[!duplicated(sumstats$chrbp) & !duplicated(sumstats$chrbp, fromLast = T),]

##check for allele code quirks (remove variants that don't have ATGC alleles and make sure those that do have capital letters)
#remove variants which have a dash for allele code, n = 9,535,795
nonstandard_alleles <- sum((!grepl("a", sumstats$Allele1) & !grepl("t", sumstats$Allele1) & !grepl("c", sumstats$Allele1) & !grepl("g", sumstats$Allele1)) | (!grepl("a", sumstats$Allele2) & !grepl("t", sumstats$Allele2) & !grepl("c", sumstats$Allele2) & !grepl("g", sumstats$Allele2)))
sumstats <- sumstats[(grepl("a", sumstats$Allele1) | grepl("t", sumstats$Allele1) | grepl("c", sumstats$Allele1) | grepl("g", sumstats$Allele1)) & (grepl("a", sumstats$Allele2) | grepl("t", sumstats$Allele2) | grepl("c", sumstats$Allele2) | grepl("g", sumstats$Allele2)),]

#make allele codes upper case
sumstats$Allele1 <- gsub("a", "A", sumstats$Allele1)
sumstats$Allele1 <- gsub("t", "T", sumstats$Allele1)
sumstats$Allele1 <- gsub("c", "C", sumstats$Allele1)
sumstats$Allele1 <- gsub("g", "G", sumstats$Allele1)
sumstats$Allele2 <- gsub("a", "A", sumstats$Allele2)
sumstats$Allele2 <- gsub("t", "T", sumstats$Allele2)
sumstats$Allele2 <- gsub("c", "C", sumstats$Allele2)
sumstats$Allele2 <- gsub("g", "G", sumstats$Allele2)

##update variant ids using the HRC SNPs file
#how about using the dbSNP file?
dbsnp37 <- fread("/home/clarklc/Programs/dbsnp37_rsonly.txt", header = F)
dbsnp37$chrbp <- paste0(dbsnp37$V1, ":", dbsnp37$V4)
dbsnp37[ ,c("V1","V3", "V4") := NULL]
names(dbsnp37)[names(dbsnp37)=="V2"] <- "SNP"
sumstats <- merge(sumstats, dbsnp37, by = "chrbp", all.x = T)

#leaves 516,808 variants without rs number

##not doing this step since all our data will now be on build 38 and 
#we don't want to mix up variants between builds
##replace the chrbp for the variants which are missing rs number
#sumstats$SNP[is.na(sumstats$SNP)] <- sumstats$chrbp[is.na(sumstats$SNP)]
sumstats <- sumstats[!is.na(sumstats$SNP),]#n variants = 9022044

##subset columns to just SNP, A1, BETA, P
sumstats <- sumstats[,c("SNP", "Allele1", "Effect", "P-value")]
names(sumstats) <- c("SNP", "A1", "BETA", "P")

##write out
write.table(sumstats, "/nfs/mahone1/PRS/testing/Kunkle_noPAC_summary_stats_cleaned.txt", row.names = F, sep = " ", quote = F)
