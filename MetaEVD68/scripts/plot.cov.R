inargs <-commandArgs(trailingOnly = TRUE)
dir <- inargs[1]
f <- inargs[2]

w=15  #image width
h=5   #image height

#Initialize plot
png(filename= paste(dir, "/basecov.png", sep=''), 
    type="cairo",
    units="in", 
    width=w, 
    height=h, 
    pointsize=12, 
    res=96)
par(pin=c(w,h), mar= c(3,5,0,1), oma = c(3, 1, 7 , 1))  #mar(bottom, left, top, right)

#nameinfo: 1=date, 2=flowcell, 3=sample, 4=index; nameinfo[[1]][X]
#filedirs <- strsplit(f, '/')
#nameinfo <- strsplit(filedirs[[1]][length(filedirs[[1]])], '_')
covData <- read.csv(f, header=TRUE, sep = '\t')[,3]

plot(covData, type= 'l', xlab = 'position', ylab = 'depth', col = 'blue')
title('Read depth per position', outer=TRUE)

dev.off()
