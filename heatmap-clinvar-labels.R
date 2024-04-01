#!/usr/bin/Rscript
library(ggplot2)
library(ggpubr)
library(BuenColors)
library(viridis)

data=read.table("clinvar-classifications-compiled.txt",sep="\t",header=TRUE)
data$gene = factor(data$gene,levels=(unique(data$gene)))
data$ClinVarLabel = factor(data$ClinVarLabel,levels=rev(c("benign","pathogenic","VUS")))

vdata = data[,c(1,2,4)] # variant-level data
pdata = data[,c(1,2,5)] # participant-level data

g1 = ggplot(vdata,aes(x=gene,y=ClinVarLabel,fill=variant_pct))+geom_tile(color="white",linewidth=1)+ 
  scale_fill_gradientn(colors = jdb_palette("brewer_celsius")) + theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",axis.title.y = element_blank()) +geom_text(aes(label=signif(variant_pct,digits=3),color = variant_pct < 10)) + coord_equal() + scale_color_manual(guide = FALSE, values = c("black", "white"))

g2 = ggplot(pdata,aes(x=gene,y=ClinVarLabel,fill=participant_pct))+geom_tile(color="white",linewidth=1)+ 
  scale_fill_gradientn(colors = jdb_palette("brewer_celsius")) + theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",axis.title.y = element_blank()) +geom_text(aes(label=signif(participant_pct,digits=3),color = participant_pct < 10)) + coord_equal() + scale_color_manual(guide = FALSE, values = c("black", "white"))

library(cowplot)
pdf("clinvar-labels.variants-participants.heatmap.pdf",height=4,width=8,useDingbats=FALSE)
plot_grid(g1,g2,ncol=2)
dev.off()

###########
mvdata = read.table("model-vs-clinvar-classifications.txt",sep="\t",header=TRUE)
mvdata$gene = factor(mvdata$gene,levels=rev(unique(mvdata$gene)))
mvdata$ClinVarLabel = factor(mvdata$ClinVarLabel,levels=c("benign_pct","pathogenic_pct","VUS_pct"))
mvdata$mergeLabel = paste(mvdata$classifier,mvdata$ModelLabel,sep=".")
#mvdata$mergeLabel = factor(mvdata$mergeLabel,levels=c("AM.benign","AM.ambiguous","AM.pathogenic","EVE.benign","EVE.ambiguous","EVE.pathogenic","ESM1b.benign","ESM1b.pathogenic"))
mvdata$mergeLabel = factor(mvdata$mergeLabel,levels=c("AM.benign","ESM1b.benign","EVE.benign","AM.ambiguous","EVE.ambiguous","AM.pathogenic","ESM1b.pathogenic","EVE.pathogenic"))

g3 = ggplot(mvdata,aes(y=gene,x=mergeLabel,fill=value,group=ClinVarLabel))+geom_tile(color="white",linewidth=1) + 
  scale_fill_gradientn(colors = jdb_palette("brewer_violet"),na.value = 'gray77') + theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",axis.title.y = element_blank()) +geom_text(aes(label=signif(value,digits=3),color = value > 70)) + geom_text(data=subset(mvdata,is.na(value)), aes(label="NA")) + coord_equal() + scale_color_manual(guide = FALSE, values = c("black", "white"))

g3 = facet(g3,facet.by="ClinVarLabel")
pdf("model-vs-clinvar-labels.variantLevel.heatmap.pdf",height=3.5,width=12,useDingbats=FALSE)
g3
dev.off()



mvdata$mergeLabel = factor(mvdata$mergeLabel,levels=c("AM.benign","AM.ambiguous","AM.pathogenic","ESM1b.benign","ESM1b.pathogenic","EVE.benign","EVE.ambiguous","EVE.pathogenic"))

g4 = ggplot(mvdata,aes(y=gene,x=mergeLabel,fill=value,group=ClinVarLabel))+geom_tile(color="white",linewidth=1) + 
  scale_fill_gradientn(colors = jdb_palette("brewer_violet"),na.value = 'gray77') + theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",axis.title.y = element_blank()) +geom_text(aes(label=signif(value,digits=3),color = value > 70)) + geom_text(data=subset(mvdata,is.na(value)), aes(label="NA")) + coord_equal() + scale_color_manual(guide = "none", values = c("black", "white"))

g4 = facet(g4,facet.by="ClinVarLabel")
pdf("model-vs-clinvar-labels.variantLevel.heatmap.reorder.pdf",height=3.5,width=12,useDingbats=FALSE)
g4
dev.off()
