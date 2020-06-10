args <- commandArgs(TRUE)
stem <- args[1]
PRSfiles <- args[2:length(args)]

library(data.table)

#get the stems for the PRS
PRS_pvals <- gsub(".profile", "", gsub(".*Pval_", "", PRSfiles))

#read in the PRS
scores <- lapply(PRSfiles, fread)

#combine into one df
ids <- scores[[1]][,c("FID", "IID")]
scores <- cbind(ids, sapply(scores, "[[", 6))
names(scores) <- c("FID", "IID", paste(stem, PRS_pvals, sep = "_"))

#write out
write.table(scores, paste0(stem, "_PRS.txt"), col.names = T, row.names = F, sep = " ", quote = F)