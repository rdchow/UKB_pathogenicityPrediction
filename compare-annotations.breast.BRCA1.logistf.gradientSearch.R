#!/usr/bin/Rscript
### Breast Cancer analysis
# BRCA1

library(data.table)
data = data.frame(fread("../brca1.sample.variant.table.scores.txt",sep="\t"),row.names=1)
diag = data.frame(fread("../../data_participant.cancerDiagnoses.simple.txt",sep="\t"),row.names=1)

dataf = data[rownames(data) %in% rownames(diag),]
diagf = diag[rownames(diag) %in% rownames(data),]

dataf2 = dataf[match(rownames(diagf),rownames(dataf)),]
all(rownames(diagf) == rownames(dataf2))

dataf2$Cancer = diagf$breast
dataf2$Sex = diagf$Sex

# binarize age by median
dataf2$AgeEnrolled = as.numeric(diagf$AgeEnrolled >= median(diagf$AgeEnrolled))

# calculate % of female participants with cancer
testdf = dataf2[ dataf2$Sex == "Female",]
table(testdf$Cancer)

# subset to females with a missense mutation and no frameshift or stop mutations
dataf3 = dataf2[dataf2$BRCA1_missense == 1 & dataf2$BRCA1_fs_stop == 0 & dataf2$Sex == "Female",]

#compare all benign to all pathogenic
library(logistf)
# include age as covariate

# ClinVar reference
model_cv = logistf(Cancer ~ ClinVar + AgeEnrolled, data=dataf3[dataf3$ClinVar %in% c("benign","pathogenic"),])

# AlphaMissense: score range 0->1
# default threshold: 0.34 (benign <-> ambig) and 0.564 (ambig <-> pathogenic)
dataf3$origAMbin = "VUS" # default status
dataf3[dataf3$AlphaMissense > 0.564,"origAMbin"] = "pathogenic"
dataf3[dataf3$AlphaMissense < 0.34,"origAMbin"] = "benign"

model_am = logistf(Cancer ~ origAMbin + AgeEnrolled, data=dataf3[dataf3$origAMbin %in% c("benign","pathogenic"),])

# gradient search for optimal thresholds
# let's try varying the pathogenic threshold to increase specificity of pathogenic calls
am_thresholds = seq(0.564:1,by=0.01)
am_params = matrix(nrow=44,ncol=8)

for (i in 1:length(am_thresholds)){
    mydata = dataf3
    mydata$AMbin = "VUS" # default status
    mydata[mydata$AlphaMissense > am_thresholds[i],"AMbin"] = "pathogenic"
    mydata[mydata$AlphaMissense < 0.34,"AMbin"] = "benign"

    confmat = as.data.frame.matrix(table(mydata$Cancer,mydata$AMbin))
    ppv = confmat[2,2] / (confmat[1,2] + confmat[2,2]) * 100

    if (nrow(mydata[mydata$AMbin == "pathogenic",]) >= 2) {
        model_am = logistf(Cancer ~ AMbin + AgeEnrolled, data=mydata[mydata$AMbin %in% c("benign","pathogenic"),])
        am_params[i,1] = am_thresholds[i]
        am_params[i,2] = as.numeric(model_am$coef[2])
        am_params[i,3] = as.numeric(confint(model_am)[2,1])
        am_params[i,4] = as.numeric(confint(model_am)[2,2])
        am_params[i,5] = as.numeric(model_am$prob[2])
        am_params[i,6] = as.numeric(table(mydata$AMbin)["VUS"])/nrow(mydata)*100
        am_params[i,7] = as.numeric(table(mydata$AMbin)["pathogenic"])/nrow(mydata)*100
        am_params[i,8] = as.numeric(ppv)
    }
    else {
        am_params[i,1] = am_thresholds[i]
        am_params[i,2:8] = rep(NA,7)
    }
}
am_params = as.data.frame(am_params)
colnames(am_params) = c("pathog_threshold","beta_coef","beta_coef.lowCI","beta_coef.hiCI","pvalue","pctVUS","pctPathogenic","PPV")
write.table(am_params,"brca1-AM-paramSearch.breast.txt",sep="\t",row.names=FALSE)


#### plot the different models based on varying thresholds
library(ggplot2)
library(ggpubr)
library(rcartocolor)
am_params = read.table("brca1-AM-paramSearch.breast.txt",sep="\t",header=TRUE)
am_params = am_params[!is.na(am_params$pvalue),]
library(dplyr)
am_params = am_params %>% mutate_at(c("pathog_threshold","beta_coef","beta_coef.lowCI","beta_coef.hiCI","pvalue","pctVUS","pctPathogenic","PPV"), as.numeric)
am_params$nlp = -log10(am_params$pvalue)
am_params[am_params$nlp > 4,"nlp"] = 4

pdf("BRCA1.breast.AMthreshold.curve.pdf",height=5,width=7,useDingbats=FALSE)
ggplot(am_params,aes(x=pathog_threshold,y=beta_coef))+geom_errorbar(aes(ymin = beta_coef.lowCI, ymax= beta_coef.hiCI), colour = "gray44",width=0) + geom_point(aes(size=nlp,fill=PPV),pch=21,color="black") +scale_fill_carto_c(palette="SunsetDark") + geom_text(data=subset(am_params,pvalue < 0.05), aes(label="*"),size=4, hjust=0.5,vjust = 0.8) + theme_bw() + scale_size_continuous(range = c(0.5, 5),breaks=c(1,2,3,4),limits=c(0,4)) + geom_hline(yintercept = as.numeric(model_cv$coef[2]), color="dodgerblue3",linetype=2) + geom_hline(yintercept = 0, color="red3",linetype=2) + annotate("rect",ymin=confint(model_cv)[2,1], ymax=confint(model_cv)[2,2],xmin=-Inf,xmax=Inf, linetype=2, fill="dodgerblue3",alpha=0.1)
dev.off()

########################################### 
### testing forest plots #####
#########################################################################
# use a higher pathogenic threshold identified from the single classifers
dataf3$AMbin = "VUS" # default status
newcut = 0.7 ### CHANGE THIS THRESHOLD MANUALLY **********

dataf3[dataf3$AlphaMissense > newcut,"AMbin"] = "pathogenic"
dataf3[dataf3$AlphaMissense < 0.34,"AMbin"] = "benign"
table(dataf3$AMbin,dataf3$ClinVar)

model_optimAM = logistf(Cancer ~ AMbin + AgeEnrolled, data=dataf3[dataf3$AMbin %in% c("benign","pathogenic"),])

######################################
# do a subanalysis looking at patients with a ClinVar VUS, then comparing AlphaMissense benign vs pathogenic
dataf4 = dataf3[dataf3$ClinVar == "VUS",]
model_cv.vus_am = logistf(Cancer ~ AMbin + AgeEnrolled, data=dataf4[dataf4$AMbin %in% c("benign","pathogenic"),])

# do another subanalysis taking ClinVar labels, then supplementing all VUS's with AlphaMissense labels
dataf3$composite_am = dataf3$ClinVar

library(tidyverse)
dataf3 = dataf3 %>% mutate(composite_am = ifelse(ClinVar == "VUS",AMbin,ClinVar)) 

model_compos.am = logistf(Cancer ~ composite_am + AgeEnrolled, data=dataf3[dataf3$composite_am %in% c("benign","pathogenic"),])


######################################
# Compile results
myres = matrix(nrow=5,ncol=10)
myres[1,] = c("ClinVar",as.numeric(model_cv$coef[2]),confint(model_cv)[2,1],confint(model_cv)[2,2],as.numeric(model_cv$prob[2]), as.numeric(table(dataf3$ClinVar)[3]/nrow(dataf3)*100) , as.numeric(model_cv$coef[3]),confint(model_cv)[3,1],confint(model_cv)[3,2],as.numeric(model_cv$prob[3]))

myres[2,] = c("AM",as.numeric(model_am$coef[2]),confint(model_am)[2,1],confint(model_am)[2,2],as.numeric(model_am$prob[2]), as.numeric(table(dataf3$origAMbin)[3]/nrow(dataf3)*100), as.numeric(model_am$coef[3]),confint(model_am)[3,1],confint(model_am)[3,2],as.numeric(model_am$prob[3]))

# ==========

myres[3,] = c("optimAM",as.numeric(model_optimAM$coef[2]),confint(model_optimAM)[2,1],confint(model_optimAM)[2,2],as.numeric(model_optimAM$prob[2]), as.numeric(table(dataf3$AMbin)[3]/nrow(dataf3)*100), as.numeric(model_optimAM$coef[3]),confint(model_optimAM)[3,1],confint(model_optimAM)[3,2],as.numeric(model_optimAM$prob[3]))

myres[4,] = c("Composite_optimAM",as.numeric(model_compos.am$coef[2]),confint(model_compos.am)[2,1],confint(model_compos.am)[2,2],as.numeric(model_compos.am$prob[2]), as.numeric(table(dataf3$composite_am)[3]/nrow(dataf3)*100) , as.numeric(model_compos.am$coef[3]),confint(model_compos.am)[3,1],confint(model_compos.am)[3,2],as.numeric(model_compos.am$prob[3]))

myres[5,] = c("CV.VUS_optimAM",as.numeric(model_cv.vus_am$coef[2]),confint(model_cv.vus_am)[2,1],confint(model_cv.vus_am)[2,2],as.numeric(model_cv.vus_am$prob[2]), as.numeric(table(dataf4$AlphaMissense)[3]/nrow(dataf4)*100) , as.numeric(model_cv.vus_am$coef[3]),confint(model_cv.vus_am)[3,1],confint(model_cv.vus_am)[3,2],as.numeric(model_cv.vus_am$prob[3]))

# ========== 
myres = as.data.frame(myres)
colnames(myres) = c("category","classifier_beta","classifier_beta.lowCI","classifier_beta.hiCI","classifier_p.value","pctVUS", "age_beta","age_beta.lowCI","age_beta.hiCI","age_p.value")

library(dplyr)
myres = myres %>% mutate_at(c("classifier_beta","classifier_beta.lowCI","classifier_beta.hiCI","classifier_p.value","pctVUS", "age_beta","age_beta.lowCI","age_beta.hiCI","age_p.value"), as.numeric)

write.table(myres,"BRCA1.female.BreastCancer.pathogenic-vs-benign.logistf.optimAM.txt",sep="\t",row.names=FALSE)

# visualize odds ratios
library(ggplot2)
myres$category = factor(myres$category,levels=rev(myres$category))

library(ggpubr)
g1 = ggplot(myres,aes(x=classifier_beta,y=category))+geom_errorbar(aes(xmin = classifier_beta.lowCI, xmax= classifier_beta.hiCI), colour = "black",width=0.2) + geom_point(size=3,color="black") + theme_bw() + xlim(-max(abs(c(myres$classifier_beta.hiCI,myres$classifier_beta.lowCI))),max(abs(c(myres$classifier_beta.hiCI,myres$classifier_beta.lowCI)))) + geom_vline(xintercept = 0, color="red3",linetype=2) + geom_vline(xintercept = myres[myres$category == "ClinVar","classifier_beta"], color="dodgerblue2",linetype=2)+ geom_text(aes(x=1,y=category,label=paste("p = ",signif(classifier_p.value,3),sep="")),nudge_y=0.25)+ geom_text(aes(x=-1,y=category,label=paste(signif(pctVUS,3),"% VUS",sep="")),nudge_y=0.25) + theme_bw() + theme(axis.line = element_line(color='black'),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())


pdf("BRCA1.female.BreastCancer.pathogenic-vs-benign.CV-vs-AM.logistf.optimAM.pdf",height=3,width=5,useDingbats=FALSE)
g1
dev.off()
