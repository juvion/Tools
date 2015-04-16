setwd("/Users/ju/Works/Sup4ochre/analysis/work19.2")
data1 <- read.table("sgle_dble_TS_delta-vs-ddG_flex37-28.in", head=T)
data2 <- read.table("sgle_dble_TS_wt-vs-ddG_flex37-28.in", head=T)
data1.1 <- subset(data1, Mutation_Count == 1)
data2.1 <- subset(data2, Mutation_Count == 1)

data1["TS_delta_norm_score"] <- (data1$delta2_reflow_normalized_GFPseq - data1$delta2_37deg_reflow_normalized_GFPseq)/data1$delta2_reflow_normalized_GFPseq 
data1.1["TS_delta_norm_score"] <- (data1.1$delta2_reflow_normalized_GFPseq - data1.1$delta2_37deg_reflow_normalized_GFPseq)/data1.1$delta2_reflow_normalized_GFPseq 
data2["TS_wt_norm_score"] <- (data2$WT2_reflow_normalized_GFPseq - data2$X37deg2_reflow_normalized_GFPseq)/data2$WT2_reflow_normalized_GFPseq
data2.1["TS_wt_norm_score"] <- (data2.1$WT2_reflow_normalized_GFPseq - data2.1$X37deg2_reflow_normalized_GFPseq)/data2.1$WT2_reflow_normalized_GFPseq

library(ggplot2)

#3D plot: dGflex_37C, TS_delta_score, delta2_reflow_normalized_GFPseq
#single_double, TS_delta
# qplot(dGflex_37C, TS_delta_score, data=data1, color = delta2_reflow_normalized_GFPseq, alpha = I(0.8)) + scale_colour_gradient(low="black", high="red") + geom_point(size=2)
# qplot(dGflex_37C, TS_wt2_score, data=data2, color = WT2_reflow_normalized_GFPseq, alpha = I(0.8)) + scale_colour_gradient(low="black", high="red") + geom_point(size=2)

# qplot(dGflex_37C, TS_delta_score, data=data1.1, color = delta2_reflow_normalized_GFPseq, alpha = I(0.8)) + scale_colour_gradient(low="black", high="red") + geom_point(size=2)
# qplot(dGflex_37C, TS_wt2_score, data=data2.1, color = WT2_reflow_normalized_GFPseq, alpha = I(0.8)) + scale_colour_gradient(low="black", high="red") + geom_point(size=2)

#plot with the normalized TS score 
pdf("dot_3D_Plot_GFP_delta_28C vs (dGflex_37C-n-TS_delta_norm).pdf", height=6, width=9)
qplot(dGflex_37C, TS_delta_norm_score, data=data1, color = delta2_reflow_normalized_GFPseq, alpha = I(0.8)) + scale_colour_gradient(low="black", high="red") + geom_point(size=2) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
dev.off()


pdf("dot_3D_Plot_GFP_delta_28C vs (dGflex_37C-n-TS_delta_norm)_single.pdf", height=6, width=9)
qplot(dGflex_37C, TS_delta_norm_score, data=data1.1, color = delta2_reflow_normalized_GFPseq, alpha = I(0.8)) + scale_colour_gradient(low="black", high="red") + geom_point(size=2) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
dev.off()


pdf("dot_3D_Plot_GFP_wt_28C vs (dGflex_37C-n-TS_wt_norm).pdf", height=6, width=9)
qplot(dGflex_37C, TS_wt_norm_score, data=data2, color = WT2_reflow_normalized_GFPseq, alpha = I(0.8)) + scale_colour_gradient(low="black", high="red") + geom_point(size=2) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
dev.off()


pdf("dot_3D_Plot_GFP_wt_28C vs (dGflex_37C-n-TS_wt_norm)_single.pdf", height=6, width=9)
qplot(dGflex_37C, TS_wt_norm_score, data=data2.1, color = WT2_reflow_normalized_GFPseq, alpha = I(0.8)) + scale_colour_gradient(low="black", high="red") + geom_point(size=2) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
dev.off()


#plot with name note
pdf("dot_3D_Plot_GFP_delta_28C vs (dGflex_37C-n-TS_delta_norm)_note.pdf", height=6, width=9)
qplot(dGflex_37C, TS_delta_norm_score, data=data1, color = delta2_reflow_normalized_GFPseq, alpha = I(0.8), label=ID) + scale_colour_gradient(low="black", high="red") + geom_point(size=2) + geom_text(size=3) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
dev.off()


pdf("dot_3D_Plot_GFP_delta_28C vs (dGflex_37C-n-TS_delta_norm)_single_note.pdf", height=6, width=9)
qplot(dGflex_37C, TS_delta_norm_score, data=data1.1, color = delta2_reflow_normalized_GFPseq, alpha = I(0.8), label=ID) + scale_colour_gradient(low="black", high="red") + geom_point(size=2) + geom_text(size=3) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
dev.off()

pdf("dot_3D_Plot_GFP_wt_28C vs (dGflex_37C-n-TS_wt_norm)_note.pdf", height=6, width=9)
qplot(dGflex_37C, TS_wt_norm_score, data=data2, color = WT2_reflow_normalized_GFPseq, alpha = I(0.8), label=ID) + scale_colour_gradient(low="black", high="red") + geom_point(size=2)+ geom_text(size=3) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
dev.off()

pdf("dot_3D_Plot_GFP_wt_28C vs (dGflex_37C-n-TS_wt_norm)_single_note.pdf", height=6, width=9)
qplot(dGflex_37C, TS_wt_norm_score, data=data2.1, color = WT2_reflow_normalized_GFPseq, alpha = I(0.8), label=ID) + scale_colour_gradient(low="black", high="red") + geom_point(size=2)+ geom_text(size=3) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
dev.off()


#ROC curve
#add a factor column to dataframe based on a conditional statement from another column
#create a new column isTS_delta based on TS_delta_score
data1$isTS <- ifelse(data1$TS_delta_score < 2, 0, 1)
pdf("ROC_Plot_dGflex_37C-n-TS_delta.pdf", height=6, width=6)
plot.roc(data1$isTS, data1$dGflex_37C, percent=F, thesholds="best", print.thres="best", 
         legacy.axes=T, xlab='false positive rate', ylab='true positive rate', main = "ROC_Plot_dGflex_37C-TS_delta")
dev.off()

data2$isTS <- ifelse(data2$TS_wt2_score < 2, 0, 1)
pdf("ROC_Plot_dGflex_37C-n-TS_wt.pdf", height=6, width=6)
plot.roc(data2$isTS, data2$dGflex_37C, percent=F, thesholds="best", print.thres="best", 
         legacy.axes=T, xlab='false positive rate', ylab='true positive rate', main = "ROC_Plot_dGflex_37C-TS_wt")
dev.off()

#ROC curve for normalized TS score
data1$isTS <- ifelse(data1$TS_delta_norm_score < 0.5, 0, 1)
pdf("ROC_Plot_dGflex_37C-n-TS_delta_norm.pdf", height=6, width=6)
plot.roc(data1$isTS, data1$dGflex_37C, percent=F, thesholds="best", print.thres="best", 
         legacy.axes=T, xlab='false positive rate', ylab='true positive rate', main = "ROC_Plot_dGflex_37C-TS_delta_norm")
dev.off()

data2$isTS <- ifelse(data2$TS_wt_norm_score < 0.5, 0, 1)
pdf("ROC_Plot_dGflex_37C-n-TS_wt_norm.pdf", height=6, width=6)
plot.roc(data2$isTS, data2$dGflex_37C, percent=F, thesholds="best", print.thres="best", 
         legacy.axes=T, xlab='false positive rate', ylab='true positive rate', main = "ROC_Plot_dGflex_37C-TS_wt_norm")
dev.off()