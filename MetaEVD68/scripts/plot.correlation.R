inargs <-commandArgs(trailingOnly = TRUE)
p_name <- inargs[1]
infile1 <- inargs[2]
infile2 <- inargs[3]

#nameinfo: 1=date, 2=flowcell, 3=sample, 4=index; nameinfo[[1]][X]
#filedirs1 <- strsplit(infile1, '/')
#nameinfo1 <- strsplit(filedirs1[[1]][length(filedirs1[[1]])], '_')
#filedirs2 <- strsplit(infile2, '/')
#nameinfo2 <- strsplit(filedirs2[[1]][length(filedirs2[[1]])], '_')


###Testing####
#infile1 <- '/Users/tanjanormark/Documents/EVD68/rastaAnalysis/bam-readcount/bam-readcount220318/merged.Aligned.RG.dedup.Q_A.bam-readcount.parsed.tsv'
#infile2 <- '/Users/tanjanormark/Documents/EVD68/rastaAnalysis/bam-readcount/bam-readcount220318/merged.Aligned.RG.dedup.Q_A.bam-readcount.parsed2.tsv'
#dir = '/Users/tanjanormark/Documents/EVD68/testing/analysis2'

library(ggplot2)
w = 10
h = 10

#Initialize plot
png(filename= p_name, 
    type="cairo",
    units="in", 
    width=w, 
    height=h, 
    pointsize=12, 
    res=150)
colors=c(0,0,139)/255
pcolor= rgb(colors[1],colors[2],colors[3], alpha=0.2)

#Get and parse data
fileData1 <- read.csv(infile1, header=TRUE, sep = '\t')[,9]
fileData2 <- read.csv(infile2, header=TRUE, sep = '\t')[,9]

#Plot
plot(log2(fileData1), log2(fileData2), main="Correlation of mismatch fractions", 
     xlab="log2 A", ylab="log2 B", pch=20, col=pcolor)
abline(lm(fileData1~fileData2), col="red") # regression line (y~x)
#Plot2
#dfData <- as.data.frame(list(fileData1, fileData2))
#colnames(dfData) <- c('SampleA','SampleB')
#
#ggplot(data=dfData, aes(x=SampleA, y=SampleB))#+

dev.off()

#Calculate statistics
corr <- cor.test(log2(fileData1), log2(fileData2), method='pearson')
results <- corr[c('estimate','p.value','statistic','method')]


#Plot########################



#abline(lm(fileData1~fileData2), col="red") # regression line (y~x)

#Plot2
#dfData <- as.data.frame(list(fileData1, fileData2))
#colnames(dfData) <- c('SampleA','SampleB')

#ggplot(data=dfData, aes(x=SampleA, y=SampleB))#+
#coord_trans(x='log2', y='log2') +
#labs(title='Correlation of mismatch fractions') +
#geom_point(color=pcolor)