args <- commandArgs(TRUE)
PRS_file <- args[1]

library(corrplot)
library(colorspace)
library(data.table)

prs <- fread(PRS_file)

#plot title
PRS_plot_filename <- gsub(".txt", ".pdf", PRS_file)

#Plot
M <- cor(prs[,names(prs)[!names(prs) %in% c("FID", "IID")], with = F], use="complete.obs")
col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
pdf(PRS_plot_filename)
corrplot(M, method="color", col=col(200),addCoef.col="black",tl.col="black",tl.srt=45)
dev.off()

#log
print(paste0("PRS correlation plot saved to ", PRS_plot_filename))