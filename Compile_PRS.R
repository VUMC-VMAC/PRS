args <- commandArgs(TRUE)
output_folder <- args[1]
genotypes_stem <- args[2]
output_tags <- args[3]
pvalues <- args[4]
output_stem <- args[5]
apoe_exclude <- args[6]

library(data.table)

setDTthreads(4)

#split the comma separated arguments
output_tags <- unlist(strsplit(output_tags, ","))
pvalues <- unlist(strsplit(pvalues, ","))

#deal with the APOE argument
apoe <- if(apoe_exclude == "yes"){ apoe <-  c("_noAPOE", "") } else { apoe <- "" }

### create list of possible profile (PRS) files ###
#use expand grid to get every possible combination (pasting additional values with genotypes_stem and pvalues to be able to paste the file name together)
PRS_files <- expand.grid(output_folder, paste0(genotypes_stem, "_"), output_tags, apoe, paste0("_PRS.Pval_", pvalues, ".profile"), stringsAsFactors = F)
names(PRS_files) <- c("output_folder", "genotypes_stem", "output_tags", "apoe", "pvalues")
#paste the file names together
PRS_files$filename <- apply(PRS_files, 1, paste0, collapse="")
#remove additional text from pvalues
PRS_files$pvalues <- gsub("_PRS.Pval_", "", gsub(".profile", "", PRS_files$pvalues))

#subset to only the PRS files that are actually present
PRS_files <- PRS_files[sapply(PRS_files$filename, file.exists),]

#read in the PRS
scores <- lapply(PRS_files$filename, fread)

#combine into one df
ids <- scores[[1]][,c("FID", "IID")]
scores <- cbind(ids, sapply(scores, "[[", 6))
names(scores) <- c("FID", "IID", paste0(PRS_files$output_tags, PRS_files$apoe, "_", PRS_files$pvalues))

#write out
write.table(scores, paste0(output_folder, "/", output_stem, ".txt"), col.names = T, row.names = F, sep = " ", quote = F)
print(paste("Scores saved to", paste0(output_folder, "/", output_stem, ".txt")))