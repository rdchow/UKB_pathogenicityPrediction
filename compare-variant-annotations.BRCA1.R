#!/usr/bin/Rscript
# Generates four plots:
# 1: ClinVar labels, on variant level
# 2: 3 stacked barplots comparing AlphaMissense, EVE, and Esm1b for recovering ClinVar labels (variant-level)
# 3: ClinVar labels, on participant level
# 4: 3 stacked barplots comparing AlphaMissense, EVE, and Esm1b for recovering ClinVar labels (participant-level)

library(data.table)
library(reshape2)
library(ggpubr)

# First, visualize ClinVar classifications at unique variant level
variants = read.table("brca1-mac4-variants-only.table.canonical.AEE.txt",sep="\t",header=TRUE)
variants = variants[grepl("missense_variant",variants$MutAnnotation),] # subset to missense variants only
variants$CLNSIG = gsub('\\.','not_provided',variants$CLNSIG) # set "." labels to "not_provided"
variants$ClinVar_clean = "" # consolidated ClinVar labels
variants[variants$CLNSIG %in% c("Benign","Benign/Likely_benign","Likely_benign"),"ClinVar_clean"] = "benign"
variants[variants$CLNSIG %in% c("Conflicting_classifications_of_pathogenicity","not_provided","Uncertain_significance"),"ClinVar_clean"] = "VUS"
variants[variants$CLNSIG %in% c("Pathogenic","Pathogenic/Likely_pathogenic","Likely_pathogenic"),"ClinVar_clean"] = "pathogenic"

# double check that the original ClinVar labels and consolidated ClinVar labels are as expected
table(variants$CLNSIG,variants$ClinVar_clean)

# Consolidate pathogenic & benign ClinVar calls, but keep the "VUS" category split up as defined by ClinVar
variants[variants$CLNSIG %in% c("Benign","Benign/Likely_benign","Likely_benign"),"CLNSIG"] = "benign"
variants[variants$CLNSIG %in% c("Pathogenic","Pathogenic/Likely_pathogenic"),"CLNSIG"] = "pathogenic"

brca1_vars.cv = as.data.frame(table(variants$ClinVar_clean,variants$CLNSIG))
brca1_vars.cv$Var1 = factor(brca1_vars.cv$Var1,levels=c("benign","VUS","pathogenic")) # this is ClinVar_clean
brca1_vars.cv$Var2 = factor(brca1_vars.cv$Var2,levels=c("benign","Conflicting_classifications_of_pathogenicity","not_provided","Uncertain_significance","pathogenic")) # this is CLNSIG
brca1_vars.cv= brca1_vars.cv[brca1_vars.cv$Freq !=0,] # drop empty rows

# plot #1
g1= ggbarplot(brca1_vars.cv,x="Var1",y="Freq",fill="Var2",color=NA,label=TRUE) + scale_fill_manual(values=c("#4C68BC","#626262","#8E8E8E","#DFDFDF","#B0291E"))+ rremove("legend") +scale_y_continuous(expand=c(0,0)) + expand_limits(y=max(brca1_vars.cv$Freq)+max(brca1_vars.cv$Freq)*0.3) + ylab("# of variants")

# Second, visualize the AlphaMissense, Eve, ESM1b labels on the Clinvar classifications
# get AM labels
brca1_vars = as.data.frame.matrix(table(variants$AlphaMissenseClassification,variants$ClinVar_clean))
brca1_vars$AM = rownames(brca1_vars)
brca1_vars$benign_pct = brca1_vars$benign/sum(brca1_vars$benign)*100
brca1_vars$pathogenic_pct = brca1_vars$pathogenic/sum(brca1_vars$pathogenic)*100
brca1_vars$VUS_pct = brca1_vars$VUS/sum(brca1_vars$VUS)*100
meltdata1 = melt(brca1_vars[,-c(1:3)])
meltdata1$classifier = "AM"

# get EVE labels
variants[variants$EVEClassification == "","EVEClassification"] = "Uncertain"
brca1_vars = as.data.frame.matrix(table(variants$EVEClassification,variants$ClinVar_clean))
brca1_vars$AM = rownames(brca1_vars)
brca1_vars$benign_pct = brca1_vars$benign/sum(brca1_vars$benign)*100
brca1_vars$pathogenic_pct = brca1_vars$pathogenic/sum(brca1_vars$pathogenic)*100
brca1_vars$VUS_pct = brca1_vars$VUS/sum(brca1_vars$VUS)*100
meltdata2 = melt(brca1_vars[,-c(1:3)])
meltdata2$classifier = "EVE"

# get ESM1b labels
brca1_vars = as.data.frame.matrix(table(variants$ESM1bClassification,variants$ClinVar_clean))
brca1_vars$AM = rownames(brca1_vars)
brca1_vars$benign_pct = brca1_vars$benign/sum(brca1_vars$benign)*100
brca1_vars$pathogenic_pct = brca1_vars$pathogenic/sum(brca1_vars$pathogenic)*100
brca1_vars$VUS_pct = brca1_vars$VUS/sum(brca1_vars$VUS)*100
meltdata3 = melt(brca1_vars[,-c(1:3)])
meltdata3$classifier = "ESM1b"

mergemelt = as.data.frame(rbind(meltdata1,meltdata2,meltdata3))
mergemelt$AM = as.character(mergemelt$AM)

# standardize labels across classifiers
mergemelt[mergemelt$AM %in% c("ambiguous","Uncertain"),"AM"] = "ambiguous" 
mergemelt[mergemelt$AM %in% c("benign","Benign"),"AM"] = "benign"
mergemelt[mergemelt$AM %in% c("Pathogenic","pathogenic"),"AM"] = "pathogenic"

mergemelt$variable = factor(mergemelt$variable,levels=c("benign_pct","VUS_pct","pathogenic_pct"))
mergemelt$AM = factor(mergemelt$AM,levels=c("benign","ambiguous","pathogenic"))

# plot #2
g2 = ggbarplot(mergemelt,x="classifier",y="value",fill="AM",color=NA,facet.by="variable",label=signif(mergemelt$value,digits=3),lab.pos="in") + scale_fill_manual(values=c("#00BDC4","#9CCB86","#F2855D")) + theme(legend.position="right") +scale_y_continuous(expand=c(0,0)) + expand_limits(y=105) + ylab("% of variants")

# ==========================
# Third, visualize ClinVar classifications at the *patient* level
data = read.table("brca1.sample.variant.table.txt",sep="\t",header=TRUE,row.names=1) # 
diag = read.table("../data_participant.cancerDiagnoses.simple.txt",sep="\t",header=TRUE,row.names=1)
data$AlphaMissense = gsub("VUS","ambiguous",data$AlphaMissense)
data$EVE = gsub("\\.","VUS",data$EVE)

dataf = data[rownames(data) %in% rownames(diag),]
diagf = diag[rownames(diag) %in% rownames(data),]
dataf2 = dataf[match(rownames(diagf),rownames(dataf)),]
all(rownames(diagf) == rownames(dataf2))

dataf3 = dataf2[dataf2$BRCA1_missense == 1,]
brca1_confusion = as.data.frame.matrix(table(dataf3$ClinVar,dataf3$AlphaMissense))
brca1_cv = as.data.frame(rowSums(brca1_confusion))
brca1_cv$category = rownames(brca1_cv)
brca1_cv$category = factor(brca1_cv$category,levels = c("benign","VUS","pathogenic"))
colnames(brca1_cv) = c("counts","category")
brca1_cv$pct = brca1_cv$counts/sum(brca1_cv$counts)*100

# plot #3
# plot the barplot with the breakdown of ClinVar classifications for each patient
g3 = ggbarplot(brca1_cv,x="category",y="pct",fill="category",color=NA,label=signif(brca1_cv$pct,digits=3)) + ylab("% of patients") + scale_fill_manual(values=c("#4C68BC","#8E8E8E","#B0291E")) + rremove("legend")+scale_y_continuous(expand=c(0,0))  + expand_limits(y=max(brca1_cv$pct)*1.1)

# Fourth, visualize AM vs EVE vs ESM1b classifications on ClinVar labels
# AlphaMissense
brca1_confusion = as.data.frame.matrix(table(dataf3$ClinVar,dataf3$AlphaMissense))
brca1_confusion_pct = brca1_confusion/rowSums(brca1_confusion)*100
brca1_confusion_pct$ClinVar = rownames(brca1_confusion_pct)
brca1_conf_melt1 = melt(brca1_confusion_pct)
brca1_conf_melt1$classifier = "AM"

# EVE
brca1_confusion = as.data.frame.matrix(table(dataf3$ClinVar,dataf3$EVE))
brca1_confusion_pct = brca1_confusion/rowSums(brca1_confusion)*100
brca1_confusion_pct$ClinVar = rownames(brca1_confusion_pct)
brca1_conf_melt2 = melt(brca1_confusion_pct)
brca1_conf_melt2$classifier = "EVE"

# ESM1b
brca1_confusion = as.data.frame.matrix(table(dataf3$ClinVar,dataf3$ESM1b))
brca1_confusion_pct = brca1_confusion/rowSums(brca1_confusion)*100
brca1_confusion_pct$ClinVar = rownames(brca1_confusion_pct)
brca1_conf_melt3 = melt(brca1_confusion_pct)
brca1_conf_melt3$classifier = "ESM1b"

# Merge all classifications into one data frame
brca1_conf_melt = as.data.frame(rbind(brca1_conf_melt1,brca1_conf_melt2,brca1_conf_melt3))
brca1_conf_melt[brca1_conf_melt$variable %in% c("VUS"),"variable"] = "ambiguous"

brca1_conf_melt$ClinVar = factor(brca1_conf_melt$ClinVar,levels=c("benign","VUS","pathogenic"))
brca1_conf_melt$variable = factor(brca1_conf_melt$variable,levels=c("benign","ambiguous","pathogenic"))
brca1_conf_melt$value = signif(brca1_conf_melt$value,digits=3)

# plot #4
g4 = ggbarplot(brca1_conf_melt,x="classifier",y="value",fill="variable",color=NA,group="classifier",label=signif(brca1_conf_melt$value,digits=3),facet.by="ClinVar",lab.vjust=1) + ylab("% of patients") + scale_fill_manual(values=c("#00BDC4","#9CCB86","#F2855D")) + theme(legend.position="right")+scale_y_continuous(expand=c(0,0)) + expand_limits(y=105)

library(cowplot)
pdf("brca1-clinvar-variant.patientLevel.annot.pdf",height=6,width=9,useDingbats=FALSE)
plot_grid(g1,g2,g3,g4,ncol=2,align="v",axis="l",rel_widths = c(0.7,2))
dev.off()
