#!/usr/bin/Rscript
### Breast Cancer analysis
# BRCA1

library(data.table)
data = read.table("brca1.sample.variant.table.txt",sep="\t",header=TRUE,row.names=1)
diag = read.table("data_participant.cancerDiagnoses.simple.txt",sep="\t",header=TRUE,row.names=1)

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

#compare all benign vs pathogenic variant carriers
library(logistf)

# include age at enrollment as a covariate
model_cv = logistf(Cancer ~ ClinVar + AgeEnrolled, data=dataf3[dataf3$ClinVar %in% c("benign","pathogenic"),])
model_am = logistf(Cancer ~ AlphaMissense + AgeEnrolled, data=dataf3[dataf3$AlphaMissense %in% c("benign","pathogenic"),])
model_eve = logistf(Cancer ~ EVE + AgeEnrolled, data=dataf3[dataf3$EVE %in% c("benign","pathogenic"),])
model_esm = logistf(Cancer ~ ESM1b + AgeEnrolled, data=dataf3[dataf3$ESM1b %in% c("benign","pathogenic"),])


# do a subanalysis looking at patients with a ClinVar VUS, then comparing AlphaMissense benign vs pathogenic
dataf4 = dataf3[dataf3$ClinVar == "VUS",]
model_cv.vus_am = logistf(Cancer ~ AlphaMissense + AgeEnrolled, data=dataf4[dataf4$AlphaMissense %in% c("benign","pathogenic"),])
model_cv.vus_eve = logistf(Cancer ~ EVE + AgeEnrolled, data=dataf4[dataf4$EVE %in% c("benign","pathogenic"),])
model_cv.vus_esm = logistf(Cancer ~ ESM1b + AgeEnrolled, data=dataf4[dataf4$ESM1b %in% c("benign","pathogenic"),])

# do another subanalysis taking ClinVar labels, then supplementing all VUS's with AlphaMissense labels
# eg, composite classifiers
dataf3$composite_am = dataf3$BRCA1_ClinVar
dataf3$composite_eve = dataf3$BRCA1_ClinVar
dataf3$composite_esm = dataf3$BRCA1_ClinVar

library(tidyverse)
dataf3 = dataf3 %>% mutate(composite_am = ifelse(ClinVar == "VUS",AlphaMissense,ClinVar)) 
dataf3 = dataf3 %>% mutate(composite_eve = ifelse(ClinVar == "VUS",EVE,ClinVar)) 
dataf3 = dataf3 %>% mutate(composite_esm = ifelse(ClinVar == "VUS",ESM1b,ClinVar)) 

model_compos.am = logistf(Cancer ~ composite_am + AgeEnrolled, data=dataf3[dataf3$composite_am %in% c("benign","pathogenic"),])
model_compos.eve = logistf(Cancer ~ composite_eve + AgeEnrolled, data=dataf3[dataf3$composite_eve %in% c("benign","pathogenic"),])
model_compos.esm = logistf(Cancer ~ composite_esm + AgeEnrolled, data=dataf3[dataf3$composite_esm %in% c("benign","pathogenic"),])
model_cv.ben.vus = logistf(Cancer ~ ClinVar + AgeEnrolled, data=dataf3[dataf3$ClinVar %in% c("benign","VUS"),])


######################################
# Compile results
myres = matrix(nrow=10,ncol=10)
myres[1,] = c("ClinVar",as.numeric(model_cv$coef[2]),confint(model_cv)[2,1],confint(model_cv)[2,2],as.numeric(model_cv$prob[2]), as.numeric(table(dataf3$ClinVar)[3]/nrow(dataf3)*100) , as.numeric(model_cv$coef[3]),confint(model_cv)[3,1],confint(model_cv)[3,2],as.numeric(model_cv$prob[3]))

myres[2,] = c("AM",as.numeric(model_am$coef[2]),confint(model_am)[2,1],confint(model_am)[2,2],as.numeric(model_am$prob[2]), as.numeric(table(dataf3$AlphaMissense)[3]/nrow(dataf3)*100) , as.numeric(model_am$coef[3]),confint(model_am)[3,1],confint(model_am)[3,2],as.numeric(model_am$prob[3]))

myres[3,] = c("EVE",as.numeric(model_eve$coef[2]),confint(model_eve)[2,1],confint(model_eve)[2,2],as.numeric(model_eve$prob[2]), as.numeric(table(dataf3$EVE)[3]/nrow(dataf3)*100) , as.numeric(model_eve$coef[3]),confint(model_eve)[3,1],confint(model_eve)[3,2],as.numeric(model_eve$prob[3]))

myres[4,] = c("ESM1b",as.numeric(model_esm$coef[2]),confint(model_esm)[2,1],confint(model_esm)[2,2],as.numeric(model_esm$prob[2]), 0 , as.numeric(model_esm$coef[3]),confint(model_esm)[3,1],confint(model_esm)[3,2],as.numeric(model_esm$prob[3]))

# ==========
myres[5,] = c("Composite_AM",as.numeric(model_compos.am$coef[2]),confint(model_compos.am)[2,1],confint(model_compos.am)[2,2],as.numeric(model_compos.am$prob[2]), as.numeric(table(dataf3$composite_am)[3]/nrow(dataf3)*100) , as.numeric(model_compos.am$coef[3]),confint(model_compos.am)[3,1],confint(model_compos.am)[3,2],as.numeric(model_compos.am$prob[3]))

myres[6,] = c("Composite_EVE",as.numeric(model_compos.eve$coef[2]),confint(model_compos.eve)[2,1],confint(model_compos.eve)[2,2],as.numeric(model_compos.eve$prob[2]), as.numeric(table(dataf3$composite_eve)[3]/nrow(dataf3)*100) , as.numeric(model_compos.eve$coef[3]),confint(model_compos.eve)[3,1],confint(model_compos.eve)[3,2],as.numeric(model_compos.eve$prob[3]))

myres[7,] = c("Composite_ESM1b",as.numeric(model_compos.esm$coef[2]),confint(model_compos.esm)[2,1],confint(model_compos.esm)[2,2],as.numeric(model_compos.esm$prob[2]), 0 , as.numeric(model_compos.esm$coef[3]),confint(model_compos.esm)[3,1],confint(model_compos.esm)[3,2],as.numeric(model_compos.esm$prob[3]))

# ========== 
myres[8,] = c("CV.VUS_AM",as.numeric(model_cv.vus_am$coef[2]),confint(model_cv.vus_am)[2,1],confint(model_cv.vus_am)[2,2],as.numeric(model_cv.vus_am$prob[2]), as.numeric(table(dataf4$AlphaMissense)[3]/nrow(dataf4)*100) , as.numeric(model_cv.vus_am$coef[3]),confint(model_cv.vus_am)[3,1],confint(model_cv.vus_am)[3,2],as.numeric(model_cv.vus_am$prob[3]))

myres[9,] = c("CV.VUS_EVE",as.numeric(model_cv.vus_eve$coef[2]),confint(model_cv.vus_eve)[2,1],confint(model_cv.vus_eve)[2,2],as.numeric(model_cv.vus_eve$prob[2]), as.numeric(table(dataf4$EVE)[3]/nrow(dataf4)*100) , as.numeric(model_cv.vus_eve$coef[3]),confint(model_cv.vus_eve)[3,1],confint(model_cv.vus_eve)[3,2],as.numeric(model_cv.vus_eve$prob[3]))

myres[10,] = c("CV.VUS_ESM1b",as.numeric(model_cv.vus_esm$coef[2]),confint(model_cv.vus_esm)[2,1],confint(model_cv.vus_esm)[2,2],as.numeric(model_cv.vus_esm$prob[2]), 0 , as.numeric(model_cv.vus_esm$coef[3]),confint(model_cv.vus_esm)[3,1],confint(model_cv.vus_esm)[3,2],as.numeric(model_cv.vus_esm$prob[3]))

# ========== 
# save regression results
myres = as.data.frame(myres)
colnames(myres) = c("category","classifier_beta","classifier_beta.lowCI","classifier_beta.hiCI","classifier_p.value","pctVUS", "age_beta","age_beta.lowCI","age_beta.hiCI","age_p.value")

library(dplyr)
myres = myres %>% mutate_at(c("classifier_beta","classifier_beta.lowCI","classifier_beta.hiCI","classifier_p.value","pctVUS", "age_beta","age_beta.lowCI","age_beta.hiCI","age_p.value"), as.numeric)

write.table(myres,"BRCA1.female.BreastCancer.pathogenic-vs-benign.logistf.txt",sep="\t",row.names=FALSE)

# visualize regression results as forest plots
library(ggplot2)
myres$category = factor(myres$category,levels=rev(myres$category))
myres$type = c(rep("standard",4),rep("composite",3),rep("CV.VUS",3))
myres$type = factor(myres$type,levels=c("standard","composite","CV.VUS"))

library(ggpubr)
myres1 = myres[myres$type == "standard",]
p1 = ggplot(myres1,aes(x=classifier_beta,y=category))+geom_errorbar(aes(xmin = classifier_beta.lowCI, xmax= classifier_beta.hiCI), colour = "black",width=0.2) + geom_point(size=3,color="black") + theme_bw() + xlim(-max(myres1$classifier_beta.hiCI),max(myres1$classifier_beta.hiCI)) + geom_vline(xintercept = 0, color="red3",linetype=2) + geom_text(aes(x=1,y=category,label=paste("p = ",signif(classifier_p.value,3),sep="")),nudge_y=0.25)+ geom_text(aes(x=-1,y=category,label=paste(signif(pctVUS,3),"% VUS",sep="")),nudge_y=0.25) + theme_bw() + theme(axis.line = element_line(color='black'),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + facet_wrap(~type, scales = 'free', nrow = 1, strip.position = 'top')

myres2 = myres[myres$type == "composite",]
p2 = ggplot(myres2,aes(x=classifier_beta,y=category))+geom_errorbar(aes(xmin = classifier_beta.lowCI, xmax= classifier_beta.hiCI), colour = "black",width=0.2) + geom_point(size=3,color="black") + theme_bw() + xlim(-max(myres2$classifier_beta.hiCI),max(myres2$classifier_beta.hiCI)) + geom_vline(xintercept = 0, color="red3",linetype=2) + geom_text(aes(x=1,y=category,label=paste("p = ",signif(classifier_p.value,3),sep="")),nudge_y=0.25)+ geom_text(aes(x=-1,y=category,label=paste(signif(pctVUS,3),"% VUS",sep="")),nudge_y=0.25) + theme_bw() + theme(axis.line = element_line(color='black'),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + facet_wrap(~type, scales = 'free', nrow = 1, strip.position = 'top')

myres3 = myres[myres$type == "CV.VUS",]
p3 = ggplot(myres3,aes(x=classifier_beta,y=category))+geom_errorbar(aes(xmin = classifier_beta.lowCI, xmax= classifier_beta.hiCI), colour = "black",width=0.2) + geom_point(size=3,color="black") + theme_bw() + xlim(-max(myres3$classifier_beta.hiCI),max(myres3$classifier_beta.hiCI)) + geom_vline(xintercept = 0, color="red3",linetype=2) + geom_text(aes(x=1,y=category,label=paste("p = ",signif(classifier_p.value,3),sep="")),nudge_y=0.25)+ geom_text(aes(x=-1,y=category,label=paste(signif(pctVUS,3),"% VUS",sep="")),nudge_y=0.25) + theme_bw() + theme(axis.line = element_line(color='black'),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + facet_wrap(~type, scales = 'free', nrow = 1, strip.position = 'top')

library(cowplot)
pdf("BRCA1.female.BreastCancer.pathogenic-vs-benign.logistf.pdf",height=6,width=4,useDingbats=FALSE)
plot_grid(p1,p2,p3,ncol=1,align="v")
dev.off()
