#R script to make a base frequency plot for the fraction of mismatches compared to the reference across the genome
library(zoo)
library(fields)

inargs <-commandArgs(trailingOnly = TRUE)
p_name <- inargs[1]
infiles <- inargs[-1]
numFiles <- length(infiles)

n= 6 #data downsampling every n'th position
ncolor= 50 #num of colors in color range
w= 18  #image width
h= 6*numFiles   #image height

plotname <- p_name
for (f in infiles){
  #nameinfo: 1=date, 2=flowcell, 3=sample, 4=index; nameinfo[[1]][X]
  filedirs <- strsplit(f, '/')
  nameinfo <- strsplit(filedirs[[1]][length(filedirs[[1]])], '_')
  plotname <- paste(plotname, nameinfo[[1]][3], sep= '_')
}

#Initialize plot
png(filename= paste(plotname,'png', sep='.'), 
    type='cairo',
    units='in', 
    width=w, 
    height=h, 
    pointsize=9 + 1*numFiles, 
    res=96)
colfunc <- colorRampPalette(c('white', 'blue'))
layout(matrix(c(1:numFiles), nrow=numFiles, byrow=TRUE))
par(pin=c(w,h), mar= c(3,5,0,1), oma = c(7, 1, 7 , 1), cex.lab=2, cex.main=3)  #mar(bottom, left, top, right)

for (f in infiles){
  filedirs <- strsplit(f, '/')
  nameinfo <- strsplit(filedirs[[1]][length(filedirs[[1]])], '_')
  fileData <- read.csv(f, header=TRUE, sep = '\t')
  
  #downsample the data for plotting purposes. Gets the maximum value for every n'th position
  frac <- rollapply(data.matrix(fileData[,9], rownames.force = NA), n, max, by = n)
  #Downsample x-ticks: get a list of every n'th position in genome
  x <- seq(0, length(frac)*n, by=n)
  image.plot(x, 1, frac, col=colfunc(ncolor), zlim = c(0,1), yaxt='n', xlab = '', ylab= nameinfo[[1]][3])
}
title("Mismatch fraction", outer=TRUE, xlab = 'Genomic position')

dev.off()