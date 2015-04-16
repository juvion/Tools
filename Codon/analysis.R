#read in the csv table, treat "#DIV/0!", "#NUM!" as NA.
rate <- read.csv("ICPfam_ORF_shift_consv-rate.csv", na.strings=c("#DIV/0!", "#NUM!"), head=T)
#omit the NA from column LogOdd.consv_rate._shift0
rate_shift0 <- na.omit(rate$LogOdd.consv_rate._shift0)
rate_shift1 <- na.omit(rate$LogOdd.consv_rate._shift1)
rate_shift2 <- na.omit(rate$LogOdd.consv_rate._shift2)
#AA rate
AA_rate_shift0 = 2.1613
AA_rate_shift1 = 0.5588
AA_rate_shift2 = 0.7187


x <- rate_shift0
z <- AA_rate_shift0
h<-hist(x, breaks=20, col="red", xlab="log(normalized rate)", xaxt='n',
        ylab="counts", main="Distribution of Log(normalized rate) at ORF shift0") 

#fit a normal curve
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)
#draw a line for mean of log(rate)
abline(v=mean(x), col="black", lwd=3)
#draw a line for 'GCAGCA' log(rate)
abline(v=z, col="green", lwd=2)

#generate a sequence of standard z scores.
z_scores <- seq(-4,4,1) * sd(x)
#convert the standard z scores to real values
z_marks <- z_scores + mean(x)
#mark the axis
axis(1, at=z_marks,labels=round(z_marks,digits=2), col.axis="blue", las=2, cex.axis=0.9, tck=-.01)
axis(1, at=z + 0.2, labels=round(z,digits=2), col.axis="green", las=0, cex.axis=0.9, tck=0, pos=10)

##############################
x <- rate_shift1
z <- AA_rate_shift1
h<-hist(x, breaks=20, col="red", xlab="log(normalized rate)", xaxt='n',
        ylab="counts", main="Distribution of Log(normalized rate) at ORF shift1") 

#fit a normal curve
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)
#draw a line for mean of log(rate)
abline(v=mean(x), col="black", lwd=3)
#draw a line for 'GCAGCA' log(rate)
abline(v=z, col="green", lwd=2)

#generate a sequence of standard z scores.
z_scores <- seq(-4,4,1) * sd(x)
#convert the standard z scores to real values
z_marks <- z_scores + mean(x)
#mark the axis
axis(1, at=z_marks,labels=round(z_marks,digits=2), col.axis="blue", las=2, cex.axis=0.9, tck=-.01)
axis(1, at=z + 0.2, labels=round(z,digits=2), col.axis="green", las=0, cex.axis=0.9, tck=0, pos=10)


##############################
x <- rate_shift2
z <- AA_rate_shift2
h<-hist(x, breaks=20, col="red", xlab="log(normalized rate)", xaxt='n',
        ylab="counts", main="Distribution of Log(normalized rate) at ORF shift2") 

#fit a normal curve
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)
#draw a line for mean of log(rate)
abline(v=mean(x), col="black", lwd=3)
#draw a line for 'GCAGCA' log(rate)
abline(v=z, col="green", lwd=2)

#generate a sequence of standard z scores.
z_scores <- seq(-4,4,1) * sd(x)
#convert the standard z scores to real values
z_marks <- z_scores + mean(x)
#mark the axis
axis(1, at=z_marks,labels=round(z_marks,digits=2), col.axis="blue", las=2, cex.axis=0.9, tck=-.01)
axis(1, at=z + 0.2, labels=round(z,digits=2), col.axis="green", las=0, cex.axis=0.9, tck=0, pos=10)

########################Global analysis

#read in the csv table, treat "#DIV/0!", "#NUM!" as NA.
rate <- read.csv("all_cp_ORF_shift_consv-rate.csv", na.strings=c("#DIV/0!", "#NUM!"), head=T)

#omit the NA from column LogOdd.consv_rate._shift0
rate_shift0 <- na.omit(rate$LogOdd.consv_rate._shift0)
rate_shift1 <- na.omit(rate$LogOdd.consv_rate._shift1)
rate_shift2 <- na.omit(rate$LogOdd.consv_rate._shift2)

#AA rate
AA_rate_shift0 = 2.1613
AA_rate_shift1 = 0.5588
AA_rate_shift2 = 0.7187


x <- rate_shift0
#plot the histogram
h<-hist(x, breaks=20, col="red", xlab="log(normalized rate)", xaxt='n',
        ylab="counts", main="Distribution of Log(normalized rate) at ORF shift0") 
#fit a normal curve
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)
#draw a line for mean of log(rate)
abline(v=mean(x), col="black", lwd=3)
#draw a line for 'GCAGCA' log(rate)
abline(v=AA_rate_shift0, col="green", lwd=2)

#generate a sequence of standard z scores.
z_scores <- seq(-4,4,1) * sd(rate_shift0)
#convert the standard z scores to real values
z_marks <- z_scores + mean(rate_shift0)
#mark the axis
axis(1, at=z_marks,labels=round(z_marks,digits=2), col.axis="blue", las=2, cex.axis=0.9, tck=-.01)
axis(1, at=AA_rate_shift0 + 0.2, labels=round(AA_rate_shift0,digits=2), col.axis="green", las=0, cex.axis=0.9, tck=0, pos=150)

##################
x <- rate_shift1
#plot the histogram
h<-hist(x, breaks=20, col="red", xlab="log(normalized rate)", xaxt='n',
        ylab="counts", main="Distribution of Log(normalized rate) at ORF shift1") 
#fit a normal curve
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)
#draw a line for mean of log(rate)
abline(v=mean(x), col="black", lwd=3)
#draw a line for 'GCAGCA' log(rate)
abline(v=AA_rate_shift1, col="green", lwd=2)

#generate a sequence of standard z scores.
z_scores <- seq(-4,4,1) * sd(rate_shift1)
#convert the standard z scores to real values
z_marks <- z_scores + mean(rate_shift1)
#mark the axis
axis(1, at=z_marks,labels=round(z_marks,digits=2), col.axis="blue", las=2, cex.axis=0.9, tck=-.01)
axis(1, at=AA_rate_shift1 + 0.2, labels=round(AA_rate_shift1,digits=2), col.axis="green", las=0, cex.axis=0.9, tck=0, pos=150)

##################
x <- rate_shift2
#plot the histogram
h<-hist(x, breaks=20, col="red", xlab="log(normalized rate)", xaxt='n',
        ylab="counts", main="Distribution of Log(normalized rate) at ORF shift2") 
#fit a normal curve
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)
#draw a line for mean of log(rate)
abline(v=mean(x), col="black", lwd=3)
#draw a line for 'GCAGCA' log(rate)
abline(v=AA_rate_shift2, col="green", lwd=2)

#generate a sequence of standard z scores.
z_scores <- seq(-4,4,1) * sd(rate_shift2)
#convert the standard z scores to real values
z_marks <- z_scores + mean(rate_shift2)
#mark the axis
axis(1, at=z_marks,labels=round(z_marks,digits=2), col.axis="blue", las=2, cex.axis=0.9, tck=-.01)
axis(1, at=AA_rate_shift2 + 0.2, labels=round(AA_rate_shift2,digits=2), col.axis="green", las=0, cex.axis=0.9, tck=0, pos=150)