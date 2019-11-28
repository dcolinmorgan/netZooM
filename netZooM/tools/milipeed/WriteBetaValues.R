rm(list=ls());

load('/proj/regeps/regep00/studies/LTCOPD/data/methylation/finalData/betas.clean_1436199620_V13.RData');
# write.table(betas.clean, file="BetaValues.txt", quote=FALSE, col.names=NA, sep="\t");
# write.table(pDat.clean, file="PhenoData.txt", quote=FALSE, row.names=F, sep="\t");

Bcase=apply(betas.clean[, pDat.clean[,51]=="case"], 1, mean);
Bcont=apply(betas.clean[, pDat.clean[,51]=="cont"], 1, mean);
Data=t(rbind(Bcase,Bcont));

write.table(Data,  file="MeanBetaValues.txt", quote=FALSE, col.names=NA, sep="\t");

Bcase=apply(betas.clean[, pDat.clean[,51]=="case"], 1, median);
Bcont=apply(betas.clean[, pDat.clean[,51]=="cont"], 1, median);
Data=t(rbind(Bcase,Bcont));

write.table(Data,  file="MedianBetaValues.txt", quote=FALSE, col.names=NA, sep="\t");
