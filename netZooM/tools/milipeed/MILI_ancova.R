
LM<-function(DEE,gene){
  # contrasts(DEE$COPD)=contr.poly(2)
  fmla <- as.formula(paste0("DEE$",gene, "~ as.factor(DEE$age_range)+as.factor(DEE$sex)+DEE$bmi+DEE$PY+as.factor(DEE$COPD)+DEE$FEV")) #*as.factor(DEE$COPD)"))
  fit <- glm(fmla, data=DEE)
  # p.value <- coef(summary(fit))[8]
  # p-age<-coef(summary(fit))[2,4]
  AIC <- AIC(fit)
  Deviance <- deviance(fit)
  cfit <- coef(summary(fit))
  fit2 <- aov(fmla, data=DEE)
  
  # df <- data.frame(gene = gene, LMp_age = cfit[[1]]$`Pr(>|t|)`[1], stringsAsFactors = F)#,p.COPD = cfit[18],p.bmi = cfit[19],p.py = cfit[20],AIC = AIC(fit), Deviance = deviance(fit), stringsAsFactors = F)
  df <- data.frame(gene = gene, LMp_fourty = cfit[32],LMp_fifty = cfit[31],LMp_sixty = cfit[34],LMp_seventy = cfit[33],aov_age=summary(fit2)[[1]][["Pr(>F)"]][1], stringsAsFactors = F)#,p.COPD = cfit[18],p.bmi = cfit[19],p.py = cfit[20],AIC = AIC(fit), Deviance = deviance(fit), stringsAsFactors = F)
  
   # results <- rbind(results, df)
  # return(results)
  return(df)
}


knitr::opts_chunk$set(echo = T, results = "hide",message=FALSE, warning=FALSE)
library(reticulate);library(stringr);library(dplyr);library(tidyr);library(data.table);library(gtools);library(ggplot2);library(plotly);library(reshape2);library(readr);library(networkD3);library(dplyr);library(igraph)

results <- data.frame(0,0,0,0,0,0)
colnames(results)<-c('gene','LMp_fourty','LMp_fifty','LMp_sixty','LMp_seventy','aov_decade')
control <- data.table:::fread('../../control.txt')
DD<-data.frame(control$number)
DD$no<-DD$control.number
DD$control.number<-NULL
DD$age<-control$age.x
DD$FEV<-control$FEV1FVC.x
DD$sex<-control$sex.x
DD$bmi<-control$BMI.x
DD$PY<-control$packyears.x
DD$no<-mixedsort(DD$no)
DD <- DD[-c(5), ]
# DD <- DD[-c(1), ]

case <- data.table:::fread('../../case.txt')
CC<-data.frame(case$number)
CC$no<-CC$case.number
CC$case.number<-NULL
CC$age<-case$age.x
CC$FEV<-case$FEV1FVC.x
CC$sex<-case$sex.x
CC$bmi<-case$BMI.x
CC$PY<-case$packyears.x
EE<-rbind(DD,CC)
EE$COPD<-factor(t(data.frame(lapply(EE$FEV, function(x) ifelse(x > .7, 'COPD', 'NO'))))) ##if(EE$FEV>.7){1}else{0}
EE$age_range<-factor(t(data.frame(lapply(EE$age, function(x) ifelse(x < 50, 'fourties', ifelse(x < 60, 'fifties', ifelse(x < 70, 'sixties', ifelse(x < 80, 'seventies', 'eighties')))))))) ##if(EE$FEV>.7){1}else{0}
# EE$age_range<-factor(t(data.frame(lapply(EE$age_range, function(x) ifelse(x < 60, 'fifties', x))))) ##if(EE$FEV>.7){1}else{0}
# EE$age_range<-factor(t(data.frame(lapply(EE$age_range, function(x) ifelse(x < 70, 'sixties', x))))) ##if(EE$FEV>.7){1}else{0}
# EE$age_range<-factor(t(data.frame(lapply(EE$age_range, function(x) ifelse(x < 80, 'seventies', x))))) ##if(EE$FEV>.7){1}else{0}


links<-data.table:::fread('../../links.txt',header=FALSE)
# CAO <- data.table:::fread('control/control_old/avg_control_old.txt',header=FALSE,drop=1)
# CAY <- data.table:::fread('control/control_young/avg_control_young.txt',header=FALSE,drop=1)
# COO <- data.table:::fread('control/control_old/avg_control_old.txt',header=FALSE,drop=1)
# COY <- data.table:::fread('control/control_young/avg_control_young.txt',header=FALSE,drop=1)
temp1 = list.files('control/lioness_output',pattern="*.txt",full.names = TRUE)
temp2 = list.files('case/lioness_output',pattern="*.txt",full.names = TRUE)
temp<-paste(c(temp1,temp2))
j=1
temp<-mixedsort(temp)
# temp<-temp[1:length(temp)-1]
for (i in 0:5){#round(dim(links)[1]/2500)){
  # temp %>% separate(A, into = 'D', extra = 'drop', remove = FALSE) %>% select(LETTERS[1:4])
  myfiles = lapply(temp,data.table::fread,header=FALSE,skip=i*2500,nrows=2500)
  # f <- function(x, pos) subset(x)
  # myfiles=lapply(temp,read_csv_chunked, DataFrameCallback$new(f), chunk_size = 2000)
  CCC<-dplyr::bind_cols(myfiles)
  # colnames(CCC)<-paste("Cont", 1:length(CCC), sep="")
  # rownames(CCC)<-paste(links$V1[1:dim(CCC)[1]],links$V2[1:dim(CCC)[1]], sep = "_")
  DDD<-transpose(CCC)
  colnames(DDD)<-paste(links$V1[1+(i*2500):((2500+(i*2500))-1)],links$V2[1+(i*2500):((2500+(i*2500))-1)], sep = "__")
  DEE<-cbind(EE,DDD)
  
  for(gene in names(DEE)[c(9:length(DEE))]){
    results[j,]<-LM(DEE,gene)
    j<-j+1}
  data.table:::fwrite(results,"SP_ANOVA.txt",sep="\t",row.names=FALSE,append=TRUE)
  
  # }
}
data.table:::fwrite(results,"SP_ANOVA.txt",sep="\t",row.names=FALSE)


 # data.table:::fwrite(results,"PA_cont_age05_b.txt",sep="\t",row.names=FALSE,col.names = TRUE)


# BF_case<-data.table:::fread('PA_case_age05_e.txt',header=TRUE)
BF_cont<-data.table:::fread('SP_ANOVA.txt',header=TRUE)
# 
# BF_case$fdr<-p.adjust(BF_case$aov_decade,method = "fdr")
BF_cont$fdr<-p.adjust(BF_cont$aov_decade,method = "fdr",n =dim(links)[1])
data.table:::fwrite(BF_cont,"SP_anova_BF.txt",sep="\t",row.names=FALSE,col.names=TRUE)
# data.table:::fwrite(BF_case,"PA_case_age05_BF.txt",sep="\t",row.names=FALSE,col.names=TRUE)

load('/proj/regeps/regep00/studies/LTCOPD/data/methylation/finalData/betas.clean_1436199620_V13.RData');
data.table:::fwrite(betas.clean,"betas.clean_1436199620_V13.txt",sep="\t",row.names=TRUE,col.names=TRUE)
