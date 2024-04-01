#!/usr/bin/Rscript
library(ggplot2)
library(ggpubr)
library(BuenColors)
library(viridis)
library(cowplot)

data=read.table("regression-results-compiled.breast.txt",sep="\t",header=TRUE)
data$nlp = log10(data$classifier_p.value)*-1

data$type = data$category
data[data$type %in% c("ClinVar","AM","EVE","ESM1b"),"type"] = "single"
data[data$type %in% c("Composite_AM","Composite_EVE","Composite_ESM1b"),"type"] = "composite"
data[data$type %in% c("CV.VUS_AM","CV.VUS_EVE","CV.VUS_ESM1b"),"type"] = "CV.VUS"

data$type = factor(data$type,levels=c("single","composite","CV.VUS"))
data$category = factor(data$category,levels=c("ClinVar","AM","ESM1b","EVE","Composite_AM","Composite_ESM1b","Composite_EVE","CV.VUS_AM","CV.VUS_ESM1b","CV.VUS_EVE")) 
data$gene = factor(data$gene,levels=rev(c("BRCA1","BRCA2","ATM","CHEK2","PALB2")))

sdata = data[data$type == "single",]
cdata = data[data$type == "composite",]
vdata = data[data$type == "CV.VUS",]

g1 = ggplot(sdata,aes(x=category,y=gene)) + geom_point(aes(fill=classifier_beta,size=nlp),pch=21,color = "black") + scale_fill_gradientn(colors = jdb_palette("brewer_spectra"),limits=c(-0.3,2)) + geom_text(data=subset(sdata,is.na(classifier_beta)), aes(label="NA")) + geom_text(data=subset(sdata,classifier_p.value < 0.05), aes(label="*"),size=8, hjust=0.5,vjust = 0.75) + scale_size_continuous(range = c(2, 12),breaks=c(2,4,6,8),limits=c(0,8)) + theme_bw() + coord_fixed() + theme(legend.position="none")

g2 = ggplot(cdata,aes(x=category,y=gene)) + geom_point(aes(fill=classifier_beta,size=nlp),pch=21,color = "black") + scale_fill_gradientn(colors = jdb_palette("brewer_spectra"),limits=c(-0.3,2)) + geom_text(data=subset(cdata,is.na(classifier_beta)), aes(label="NA")) + geom_text(data=subset(cdata,classifier_p.value < 0.05), aes(label="*"),size=8, hjust=0.5,vjust = 0.75) + scale_size_continuous(range = c(2, 12),breaks=c(2,4,6,8),limits=c(0,8)) + theme_bw() + coord_fixed() + theme(legend.position="none")

g3 = ggplot(vdata,aes(x=category,y=gene)) + geom_point(aes(fill=classifier_beta,size=nlp),pch=21,color = "black") + scale_fill_gradientn(colors = jdb_palette("brewer_spectra"),limits=c(-0.3,2)) + geom_text(data=subset(vdata,is.na(classifier_beta)), aes(label="NA")) + geom_text(data=subset(vdata,classifier_p.value < 0.05), aes(label="*"),size=8, hjust=0.5,vjust = 0.75) + scale_size_continuous(range = c(2, 12),breaks=c(2,4,6,8),limits=c(0,8)) + theme_bw() + coord_fixed() 

#g1 = facet(g1,facet.by="type",scales="fixed")

##### ovarian
data=read.table("regression-results-compiled.ovary.txt",sep="\t",header=TRUE)
data$nlp = log10(data$classifier_p.value)*-1

data$type = data$category
data[data$type %in% c("ClinVar","AM","EVE","ESM1b"),"type"] = "single"
data[data$type %in% c("Composite_AM","Composite_EVE","Composite_ESM1b"),"type"] = "composite"
data[data$type %in% c("CV.VUS_AM","CV.VUS_EVE","CV.VUS_ESM1b"),"type"] = "CV.VUS"

data$type = factor(data$type,levels=c("single","composite","CV.VUS"))
data$category = factor(data$category,levels=c("ClinVar","AM","ESM1b","EVE","Composite_AM","Composite_ESM1b","Composite_EVE","CV.VUS_AM","CV.VUS_ESM1b","CV.VUS_EVE")) 
data$gene = factor(data$gene,levels=rev(c("BRCA1","BRCA2")))

sdata = data[data$type == "single",]
cdata = data[data$type == "composite",]
vdata = data[data$type == "CV.VUS",]

g4 = ggplot(sdata,aes(x=category,y=gene)) + geom_point(aes(fill=classifier_beta,size=nlp),pch=21,color = "black") + scale_fill_gradientn(colors = jdb_palette("brewer_spectra"),limits=c(-0.3,3)) + geom_text(data=subset(sdata,is.na(classifier_beta)), aes(label="NA")) + geom_text(data=subset(sdata,classifier_p.value < 0.05), aes(label="*"),size=8, hjust=0.5,vjust = 0.75) + scale_size_continuous(range = c(2, 12),breaks=c(2,4,6),limits=c(0,6)) + theme_bw()  + theme(legend.position="none") + coord_fixed()

g5 = ggplot(cdata,aes(x=category,y=gene)) + geom_point(aes(fill=classifier_beta,size=nlp),pch=21,color = "black") + scale_fill_gradientn(colors = jdb_palette("brewer_spectra"),limits=c(-0.3,3)) + geom_text(data=subset(cdata,is.na(classifier_beta)), aes(label="NA")) + geom_text(data=subset(cdata,classifier_p.value < 0.05), aes(label="*"),size=8, hjust=0.5,vjust = 0.75) + scale_size_continuous(range = c(2, 12),breaks=c(2,4,6),limits=c(0,6)) + theme_bw()  + theme(legend.position="none") + coord_fixed()

g6 = ggplot(vdata,aes(x=category,y=gene)) + geom_point(aes(fill=classifier_beta,size=nlp),pch=21,color = "black") + scale_fill_gradientn(colors = jdb_palette("brewer_spectra"),limits=c(-0.3,3)) + geom_text(data=subset(vdata,is.na(classifier_beta)), aes(label="NA")) + geom_text(data=subset(vdata,classifier_p.value < 0.05), aes(label="*"),size=8, hjust=0.5,vjust = 0.75) + scale_size_continuous(range = c(2, 12),breaks=c(2,4,6),limits=c(0,6)) + theme_bw() + coord_fixed()


pdf("dotplot-regression-breast.pdf",height=7,width=13,useDingbats=FALSE)
#plot_grid(g1,g2,g3,ncol=3,align="hv",widths=c(1.2,1,1),axis="t")
ggarrange(g1,g2,g3,ncol=3,align="v",widths=c(1.153,1,1))
dev.off()


pdf("dotplot-regression-ovary.pdf",height=7,width=13,useDingbats=FALSE)
#plot_grid(g4,g5,g6,ncol=3,align="v",widths=c(1.2,1,2),axis="t")
ggarrange(g4,g5,g6,ncol=3,align="v",widths=c(1.153,1,1))
dev.off()

####### now plot the % of VUSs remaining with each classifier
library(ggplot2)
library(ggpubr)
library(BuenColors)
library(viridis)
setwd("C:/Parikh/UKBB/_meta-visualization")
data=read.table("regression-results-compiled.breast.txt",sep="\t",header=TRUE)
data$type = data$category
data[data$type %in% c("ClinVar","AM","EVE","ESM1b"),"type"] = "single"
data[data$type %in% c("Composite_AM","Composite_EVE","Composite_ESM1b"),"type"] = "composite"
data[data$type %in% c("CV.VUS_AM","CV.VUS_EVE","CV.VUS_ESM1b"),"type"] = "CV.VUS"

data$type = factor(data$type,levels=c("single","composite","CV.VUS"))
data$category = factor(data$category,levels=c("ClinVar","AM","ESM1b","EVE","Composite_AM","Composite_ESM1b","Composite_EVE","CV.VUS_AM","CV.VUS_ESM1b","CV.VUS_EVE")) 
data$gene = factor(data$gene,levels=rev(c("BRCA1","BRCA2","ATM","CHEK2","PALB2")))

sdata = data[data$type == "single",]
cdata = data[data$type == "composite",]
vdata = data[data$type == "CV.VUS",]

g7 = ggplot(sdata,aes(x=category,y=gene)) + geom_tile(aes(fill=pctVUS),color = "white",linewidth=1) + scale_fill_gradientn(colors = rev(jdb_palette("calma_morado")),limits=c(0,100)) + geom_text(data=subset(sdata,is.na(pctVUS)), aes(label="NA")) + geom_text(data=subset(sdata,!is.na(pctVUS)), aes(label=signif(pctVUS,digits=3),color = pctVUS > 50)) + theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",axis.title.y = element_blank())  + coord_fixed() + scale_color_manual(guide = "none", values = c("black", "white"))

g8 = ggplot(cdata,aes(x=category,y=gene)) + geom_tile(aes(fill=pctVUS),color = "white",linewidth=1) + scale_fill_gradientn(colors = rev(jdb_palette("calma_morado")),limits=c(0,100)) + geom_text(data=subset(cdata,is.na(pctVUS)), aes(label="NA")) + geom_text(data=subset(cdata,!is.na(pctVUS)), aes(label=signif(pctVUS,digits=3),color = pctVUS > 50)) + theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",axis.title.y = element_blank())  + coord_fixed() + scale_color_manual(guide = "none", values = c("black", "white"))

g9 = ggplot(vdata,aes(x=category,y=gene)) + geom_tile(aes(fill=pctVUS),color = "white",linewidth=1) + scale_fill_gradientn(colors = rev(jdb_palette("calma_morado")),limits=c(0,100)) + geom_text(data=subset(vdata,is.na(pctVUS)), aes(label="NA")) + geom_text(data=subset(vdata,!is.na(pctVUS)), aes(label=signif(pctVUS,digits=3),color = pctVUS > 50)) + theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",axis.title.y = element_blank())  + coord_fixed() + scale_color_manual(guide = "none", values = c("black", "white"))

library(cowplot)
pdf("pctVUS.heatmap.pdf",height=7,width=13,useDingbats=FALSE)
plot_grid(g7,g8,g9,ncol=3,align="v",widths=c(1.2,1,1),axis="t")
dev.off()
