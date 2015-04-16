#load library
library(ggplot2)
library(pROC)

#setup working dir.
workdir <- "/Users/ju/Works/Sup4ochre/analysis/work21.1"
infile <- "tRNA_1_20_15.csv"
#get the file name without extension
file_prefix <- strsplit(infile, '\\.')[[1]][1]

setwd(workdir)

#read in original csv file as raw_data
raw_data <- read.csv(infile, header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)

#########################
###Filtering data
########################

#setup filtering settings
f_reads = 100
f_cell_counts = 30
f_mutation_counts = 3
f_GFP_seq = 0.056
f_isTS_delta = 0.50
f_isTS_wt = 0.65
f_isRTD = 2.0

#subsetters
sub_mutation = 1
sub_dG_flex = -22.5 
sub_TS_norm_delta = 0.50
sub_TS_norm_wt = 0.65


#a function to test the variant is bimodal or not. return True or False.
bimodal_test <- function(bin1, bin2, bin3, bin4) {
    isBimodal <- ((bin1 > bin4) & (bin4 > 0.2*bin1) & (bin3 < 0.5*bin4 | bin2 < 0.5*bin4)) | ((bin4 > bin1) & (bin1 > 0.2*bin4) & (bin3 < 0.5*bin1 | bin2 < 0.5*bin1))
    return(isBimodal)
}

#select come columns from raw_data
collist <- c("ID", "WT2.reflow.normalized.GFPseq", "delta2.reflow.normalized.GFPseq", "X37deg.4.reflow.normalized.GFPseq", 
             "delta2.37deg.reflow.normalized.GFPseq", "Mutation.Count",
             "WT2.Total.Counts", "delta2.Total.Counts", "X37deg.4.Total.Counts", "delta2.37deg.Total.Counts",  
             "WT2.Total.Cells", "delta2.Total.Cells",  "X37deg.4.Total.Cells", "delta2.37deg.Total.Cells",
             "ID2", "pos", "mut", "loc", "mature_ddG_rigid_22C", "mature_ddG_flex_22C", "mature_ED_22C", 
             "mature_ddG_rigid_28C", "mature_ddG_flex_28C", "mature_ED_28C", "mature_ddG_rigid_37C", "mature_ddG_flex_37C", 
             "mature_ED_37C", "mature_ddG_rigid_45C", "mature_ddG_flex_45C", "mature_ED_45C", "dGflex_28C", "dGflex_37C", 
             "ddG_flex37.28", "Sequence")

data <- raw_data[collist]

#filter the data: 
filtered_data <- subset(data, Mutation.Count <= f_mutation_counts & 
                            # filtered_data <- subset(data, T & 
                            WT2.Total.Counts >= f_reads & delta2.Total.Counts >= f_reads & 
                            X37deg.4.Total.Counts >= f_reads & delta2.37deg.Total.Counts >= f_reads &  
                            WT2.Total.Cells >= f_cell_counts & delta2.Total.Cells >= f_cell_counts & 
                            X37deg.4.Total.Cells >= f_cell_counts & delta2.37deg.Total.Cells >= f_cell_counts &
                            !grepl("loop", loc) & !grepl("mismatch", loc))

write.csv(filtered_data, paste(file_prefix, 'filtered_data.csv', sep = '_'), row.names=FALSE)

#calculate TS normalized scores for delta and wt, and create a binary variable isTS
filtered_data["TS_delta_norm_score"] <- (filtered_data$delta2.reflow.normalized.GFPseq - 
                                             filtered_data$delta2.37deg.reflow.normalized.GFPseq)/filtered_data$delta2.reflow.normalized.GFPseq

filtered_data$isTS_delta <- ifelse(filtered_data$TS_delta_norm_score < f_isTS_delta, 0, 1)

filtered_data["TS_wt_norm_score"] <- (filtered_data$WT2.reflow.normalized.GFPseq - 
                                          filtered_data$X37deg.4.reflow.normalized.GFPseq)/filtered_data$WT2.reflow.normalized.GFPseq

filtered_data$isTS_wt <- ifelse(filtered_data$TS_wt_norm_score < f_isTS_wt, 0, 1)


#calculate RTD normalized scores for delta and wt, and create a binary variable isTS
filtered_data["RTD_28C_score"] <- filtered_data$delta2.reflow.normalized.GFPseq / filtered_data$WT2.reflow.normalized.GFPseq

filtered_data$isRTD_28C <- ifelse(filtered_data$RTD_28C_score < f_isRTD, "notRTD", "RTD")

filtered_data["RTD_37C_score"] <- filtered_data$delta2.37deg.reflow.normalized.GFPseq / filtered_data$X37deg.4.reflow.normalized.GFPseq

filtered_data$isRTD_37C <- ifelse(filtered_data$RTD_37C_score < f_isRTD, "notRTD", "RTD")



#filter for delta batch, using cutoff of GPF_seq
filtered_data_delta <- subset(filtered_data, delta2.reflow.normalized.GFPseq >= f_GFP_seq)
filtered_data_wt <- subset(filtered_data, WT2.reflow.normalized.GFPseq >= f_GFP_seq)


#pickout subset of the filtered_data for homogeneous stems.
#and generate a stacking stems(at_stem: acceptor and t stems, cd_stem: anticodon and d stems.).
for (batch in c('delta', 'wt')) {
    data_name <- paste('filtered_data', batch, sep='_')
    data <- get(data_name)
    for (stem in c('as', 'ds', 'cs', 'ts', 'at', 'cd')) {
        if (stem == 'as') {
            stem_name = 'acceptor_stem'
            subdata <- subset(data, grepl('acceptor', loc) & 
                                  (!grepl('d_stem', loc) & !grepl('anticodon', loc) & !grepl('t_stem', loc)))
        }
        if (stem == 'ds') {
            stem_name = 'd_stem'
            subdata <- subset(data, grepl('d_stem', loc) & 
                                  (!grepl('acceptor_stem', loc) & !grepl('anticodon', loc) & !grepl('t_stem', loc)))
        }
        if (stem == 'cs') {
            stem_name = 'anticodon_stem'
            subdata <- subset(data, grepl('anticodon_stem', loc) & 
                                  (!grepl('acceptor_stem', loc) & !grepl('d_stem', loc) & !grepl('t_stem', loc)))
        }
        if (stem == 'ts') {
            stem_name = 't_stem'
            subdata <- subset(data, grepl('t_stem', loc) & 
                                  (!grepl('acceptor_stem', loc) & !grepl('d_stem', loc) & !grepl('anticodon_stem', loc)))
        }
        #stacking stems: acceptor-t stems (at) and anticodon-d stems (ct)
        if (stem == 'at') {
            stem_name = 'at_stack_stem'
            subdata <- subset(data, (grepl('t_stem', loc) | grepl('acceptor', loc)) & 
                                  (!grepl('d_stem', loc) & !grepl('anticodon_stem', loc)))
        }        
        if (stem == 'cd') {
            stem_name = 'cd_stack_stem'
            subdata <- subset(data, (grepl('d_stem', loc) | grepl('anticodon', loc)) & 
                                  (!grepl('t_stem', loc) & !grepl('acceptor', loc)))
        }           
        subdata_name <- paste('filtered_data', batch, stem, sep='_')
        assign(subdata_name, subdata)
        write.csv(subdata, paste(subdata_name, '.csv', sep = ''), row.names=FALSE)
    }
}

#pickout subset of the filtered_data containing 'acceptor stem'.
for (batch in c('delta', 'wt')) {
    data_name <- paste('filtered_data', batch, sep='_')
    data <- get(data_name)
    subdata <- subset(data, grepl('acceptor', loc))
    subdata_name <- paste('filtered_data', batch, 'incld_as', sep='_')
    assign(subdata_name, subdata)
    write.csv(subdata, paste(subdata_name, '.csv', sep = ''), row.names=FALSE)
}


#pickout subset of the filtered_data not containing 'acceptor stem'.
for (batch in c('delta', 'wt')) {
    data_name <- paste('filtered_data', batch, sep='_')
    data <- get(data_name)
    subdata <- subset(data, !grepl('acceptor', loc))
    subdata_name <- paste('filtered_data', batch, 'excld_as', sep='_')
    assign(subdata_name, subdata)
    write.csv(subdata, paste(subdata_name, '.csv', sep = ''), row.names=FALSE)
}

#generate a combined data for selected stems (as, ds, cs, ts).
filtered_data_delta_same_stem <- rbind(filtered_data_delta_as, filtered_data_delta_ds, filtered_data_delta_cs, filtered_data_delta_ts) 
filtered_data_wt_same_stem <- rbind(filtered_data_wt_as, filtered_data_wt_ds, filtered_data_wt_cs, filtered_data_wt_ts)




######################################################################################
##### Plot scatter: variables are RTD_37C_score, ED_37C, normalized TS score for 
#### only acceptor stem region
######################################################################################



#iterate batches
for (batch in c('delta', 'wt')) {
    #iterate stems
    if (batch == 'delta') y_cutoff = 0.5
    if (batch == 'wt') y_cutoff = 0.65
    for (stem in c('', 'as', 'ds', 'cs', 'ts', 'incld_as', 'excld_as', 'same_stem')) {
        if (stem == '') stem_name = ''
        if (stem == 'as') stem_name = 'acceptor_stem'
        if (stem == 'ds') stem_name = 'd_stem'
        if (stem == 'cs') stem_name = 'anticodon_stem'
        if (stem == 'ts') stem_name = 't_stem'
        if (stem == 'incld_as') stem_name = 'incld_acceptor_stem'
        if (stem == 'excld_as') stem_name = 'excld_acceptor_stem'
        if (stem == 'same_stem') stem_name = 'same_stem'
        
        #define data_name, plot_name(with note and without note)
        if (stem == '') {
            data_name <- paste('filtered_data', batch, sep='_')
            plot_file_name1 <- paste('RTD_37C_score vs (dGflex_37C-n-TS_', batch, '.pdf', sep='')
            plot_file_name2 <- paste('RTD_37C_score vs (dGflex_37C-n-TS_', batch, '_note.pdf', sep='')
        }
        else {data_name <- paste('filtered_data', batch, stem, sep='_')
              plot_file_name1 <- paste('RTD_37C_score vs (dGflex_37C-n-TS_', batch, '_', stem, '.pdf', sep='')
              plot_file_name2 <- paste('RTD_37C_score vs (dGflex_37C-n-TS_', batch, '_', stem, '_note.pdf', sep='')
        }
        
        #get() convert the string to data name.
        plot_data <- get(data_name)
        
        #print scatter plot pdfs without id note.
        pdf(plot_file_name1, width = 9, height = 6)
        #in for loop, auto printing does not work.
        print(
            ggplot(data=plot_data, aes(dGflex_37C, get(paste('TS_', batch, '_norm_score', sep='')))) + 
                geom_point(aes(color = isRTD_37C), size=2) +  
                labs(title = paste('dGflex_37C vs ', paste('TS_', batch, '_norm_score', sep=''), ' ', stem_name, sep = ''),  
                     x = 'dGflex_37C', 
                     y = paste('TS_', batch, '_norm_score', sep='')) + 
                theme(axis.text=element_text(size=12),
                      axis.title=element_text(size=14,face="bold")) + 
                xlim(c(-30, -15)) + ylim(c(-0.2, 1)) +
                geom_vline(xintercept=c(-25,-22.5), linetype="dotted", lw = 2, col = 'blue') +
                geom_hline(yintercept=y_cutoff, linetype="dotted", lw = 2, col = 'green')
        )
        dev.off()  
        
        #print scatter plot pdfs with id note.
        pdf(plot_file_name2, width = 9, height = 6)
        #in for loop, auto printing does not work.
        print(
            qplot(dGflex_37C, get(paste(paste('TS_', batch, '_norm_score', sep=''))), 
                  data=plot_data, color =isRTD_37C , alpha = I(0.8), label=ID) + 
                labs(title = paste('dGflex_37C vs ', paste('TS_', batch, '_norm_score', sep=''), ' ', stem_name,  sep = ''),  
                     x = 'dGflex_37C', 
                     y = paste('TS_', batch, '_norm_score', sep='')) + 
                geom_point(aes(color = isRTD_37C), size=2) + 
                geom_point(size=2) + 
                geom_text(size=3) + 
                theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + 
                xlim(c(-30, -15)) + ylim(c(-0.2, 1)) +
                geom_vline(xintercept=c(-25,-22.5), linetype="dotted", lw =2, col = 'blue') +
                geom_hline(yintercept=y_cutoff, linetype="dotted", lw = 2, col = 'green')
        )
        dev.off()
    }
}








#####----------------------------------------


pick_outlier <- FALSE
if (pick_outlier == TRUE) {
    #######################
    ### pick outliers #####
    #######################
    #get subset of data
    #Get subset of variants:filtered_data_delta (single, dGlfex_37C <=-22.5kcal/mol, TS_norm_score_delta >= 0.5) 
    #                       OR (dGlfex_37C > -22.5kcal/mol, TS_norm_score_delta < 0.5)
    outlier_data_delta <- subset(filtered_data_delta, (Mutation.Count <= sub_mutation & dGflex_37C <= -22.5 &
                                                           TS_delta_norm_score >= 0.5) |(dGflex_37C > -22.5 & TS_delta_norm_score < 0.5))
    write.csv(outlier_data_delta, "tRNA_7_31_14_Ju_processed_10-29-14b_outlier_data_delta.csv", row.names=FALSE)
    write.table(outlier_data_delta[,c("ID", "Sequence")],"tRNA_7_31_14_Ju_processed_10-29-14b_outlier_data_delta_seq.in", 
                row.names=F, sep="\t", quote=F)
    
    
    #Get subset of variants:filtered_data_wt (single, dGlfex_37C <=-22.5kcal/mol, TS_norm_score_delta >= 0.65) 
    #                       OR (dGlfex_37C > -22.5kcal/mol, TS_norm_score_delta < 0.65)
    
    outlier_data_wt <- subset(filtered_data_wt, (Mutation.Count <= sub_mutation & dGflex_37C <= -22.5 &
                                                     TS_wt_norm_score >= 0.65) |(dGflex_37C > -22.5 & TS_wt_norm_score < 0.65))
    write.csv(outlier_data_wt, "tRNA_7_31_14_Ju_processed_10-29-14b_outlier_data_wt.csv", row.names=FALSE)
    write.table(outlier_data_wt[,c("ID", "Sequence")],"tRNA_7_31_14_Ju_processed_10-29-14b_outlier_data_wt_seq.in", 
                row.names=F, sep="\t", quote=F)
    
}




