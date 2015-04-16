#to run the script in terminal
#Rscript woring13.3.R 
#or
#R CMD batch woring13.3.R 

library(dplyr)

setwd("/Users/ju/Works/Codon/analysis/work13.3")

#import the dataset
ICP_z_data <- read.csv("ICP_fam_CodonPairs_zscore.csv", head=T)

# Number of replicas
n = 1000

# Level of significance
alpha = 0.05

#Difference between t tests of bootstrapped dataset (subsets isICP and noICP)
ttest.bootstrap = NULL

bootstrap = 'pool'

if (bootstrap = 'pool') {
    a <- ICP_z_data
    for (i in 1:n) {
        #Sample with replacement
        a.bootstrap = sample_n(a, size = nrow(a), replace = TRUE)
    #     ttest.bootstrap[i] = mean(a.bootstrap$LogOdd.consv_rate._shift0)
        ttest.bootstrap[i] = t.test(a.bootstrap[a.bootstrap$isICP == 0, ]$LogOdd.consv_rate._shift0, a.bootstrap[a.bootstrap$isICP == 1, ]$LogOdd.consv_rate._shift0)$statistic
    }
}

if(bootstrap != 'pool' ) {
    a <- ICP_z_data[ICP_z_data$isICP == 0,]
    b <- ICP_z_data[ICP_z_data$isICP == 1,]
    # Number of replicas
    n = 1000
    
    # Level of significance
    alpha = 0.05
    
    #Difference between t tests of bootstrapped dataset (subsets isICP and noICP)
    ttest.bootstrap = NULL
    
    for (i in 1:n) {
        #Sample with replacement
        a.bootstrap = sample_n(a, size = nrow(a), replace = TRUE)
        b.bootstrap = sample_n(b, size = nrow(a), replace = TRUE)
        #     ttest.bootstrap[i] = mean(a.bootstrap$LogOdd.consv_rate._shift0)
        ttest.bootstrap[i] = t.test(a.bootstrap$LogOdd.consv_rate._shift0, b.bootstrap$LogOdd.consv_rate._shift0)$statistic
    }
}

# Confidence interval
quantile(ttest.bootstrap, c(alpha/2, 1 - alpha/2))
mean(ttest.bootstrap)

##############################################
###plot codon pairs conservation histogram####
##############################################

ggplot(ICP_z_data, aes(x=Z_LogOdd.consv_rate._shift0, fill=as.character(isICP))) + 
  geom_histogram(binwidth=.5, alpha=.5, position="identity") +
  labs(x='conservation rate Z-score', y='count', title = 'ICP Dipeptide Famliy Z-score Distribution\n ORF shift0')
ggsave(filename = 'ICP_z_score_hist_shift0.png')

ggplot(ICP_z_data, aes(x=Z_LogOdd.consv_rate._shift1, fill=as.character(isICP))) + 
  geom_histogram(binwidth=.5, alpha=.5, position="identity") +
  labs(x='conservation rate Z-score', y='count', title = 'ICP Dipeptide Famliy Z-score Distribution\n ORF shift1')
ggsave(filename = 'ICP_z_score_hist_shift1.png')


ggplot(ICP_z_data, aes(x=Z_LogOdd.consv_rate._shift2, fill=as.character(isICP))) + 
  geom_histogram(binwidth=.5, alpha=.5, position="identity") +
  labs(x='conservation rate Z-score', y='count', title = 'ICP Dipeptide Famliy Z-score Distribution\n ORF shift2')
ggsave(filename = 'ICP_z_score_hist_shift2.png')
