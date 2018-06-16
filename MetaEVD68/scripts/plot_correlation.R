inargs <-commandArgs(trailingOnly = TRUE)
p_name <- inargs[1]
infile1 <- inargs[2]
infile2 <- inargs[3]

#nameinfo: 1=date, 2=flowcell, 3=sample, 4=index; nameinfo[[1]][X]
filedirs1 <- strsplit(infile1, '/')
nameinfo1 <- strsplit(filedirs1[[1]][length(filedirs1[[1]])], '_')
filedirs2 <- strsplit(infile2, '/')
nameinfo2 <- strsplit(filedirs2[[1]][length(filedirs2[[1]])], '_')

sample1 <- nameinfo1[[1]][3]
sample2 <- nameinfo2[[1]][3]

w = 10
h = 10

#Initialize plot
png(filename= paste(paste(p_name,nameinfo1[[1]][3],nameinfo2[[1]][3], sep='_'),'png', sep='.'), 
    type="cairo",
    units="in", 
    width=w, 
    height=h, 
    pointsize=12, 
    res=150)
colors=c(0,0,139)/255
pcolor= rgb(colors[1],colors[2],colors[3], alpha=0.2)
par(cex.lab=1.5, cex.main=2, oma = c(1,1,1,1)) #mar(bottom, left, top, right)

#Get and parse data
fileData1 <- read.csv(infile1, header=TRUE, sep = '\t')[,9]
fileData2 <- read.csv(infile2, header=TRUE, sep = '\t')[,9]
df <- do.call(rbind, Map(data.frame, A=fileData1, B=fileData2))
filter_df <- subset(df , A > 0.01 | B > 0.01)

#Plot
x_data <- (filter_df[,1] + 10^(-5))
y_data <- (filter_df[,2] + 10^(-5))
plot(x_data, y_data, main="Correlation of variations (iSNVs) per genomic position", 
     xlab=paste("log10 fraction of variants", sample1), ylab=paste("log10 fraction of variants", sample2), pch=20, col=pcolor, log = 'xy')
abline(lm(log10(y_data)~log10(x_data)), col="red") # regression line (y~x)
abline(v=0.01,col="black")
abline(h=0.01,col="black")

dev.off()

#Calculate statistics
corr <- cor.test(filter_df[,1], filter_df[,2], method='pearson')
results <- corr[c('estimate','p.value','statistic','method', 'conf.int')]
print('Correlation_results:')
print(results)
print('Correlation_end')
print(sample1)
print(sum(fileData1>0.01)/length(fileData1))
print(length(fileData1))
print(sum(fileData1>0.01))
print(sample2)
print(length(fileData2))
print(sum(fileData2>0.01))
print(sum(fileData2>0.01)/length(fileData2))
