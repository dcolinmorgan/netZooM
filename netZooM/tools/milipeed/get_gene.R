library(illuminaHumanv4.db)
load("/proj/regeps/regep00/studies/LTCOPD/data/expression/finalData/expSet.clean_V14_1435701821.RData")
genes <- read.table("~/analyses/MILIPEED/genes.txt", quote="\"", comment.char="")
w<-data.frame(Gene=unlist(mget(x = genes$V1,envir = illuminaHumanv4SYMBOL)))
# is.na(w$Gene)
# exp.clean[FALSE==is.na(w$Gene),]
z<-exp.clean[FALSE==is.na(w$Gene),]
# rownames(z)<-w$Gene
# rownames(z)<-w$Gene[FALSE==is.na(w$Gene),]
rownames(z)<-w$Gene[FALSE==is.na(w$Gene)]
z<-data.frame(z)
data.table:::fwrite(z,"gene_exp_v14b.txt",sep="\t",row.names=TRUE,col.names=TRUE)
write.csv2(rownames(z),"genes2.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote = FALSE)
