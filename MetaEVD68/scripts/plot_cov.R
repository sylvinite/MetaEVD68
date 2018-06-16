inargs <-commandArgs(trailingOnly = TRUE)
f_dir <- inargs[1]
f <- inargs[2]

w=15  #image width
h=5   #image height

#Initialize plot
#nameinfo: 1=date, 2=flowcell, 3=sample, 4=index; nameinfo[[1]][X]
filedirs <- strsplit(f, '/')
filedirs_info <- strsplit(filedirs[[1]][length(filedirs[[1]])],'[.]')
nameinfo <- strsplit(filedirs[[1]][length(filedirs[[1]])], '_')
f_name = paste(f_dir, '/', filedirs_info[[1]][1], ".cov.png", sep='')
print(paste('FILENAME:',f_name,':ENDFILE', sep=''))

png(filename= f_name, 
    type="cairo",
    units="in", 
    width=w, 
    height=h, 
    pointsize=12, 
    res=96)
par(pin=c(w,h), mar= c(3,5,0,1), oma = c(3, 1, 7 , 1))  #mar(bottom, left, top, right)


covData <- read.csv(f, header=TRUE, sep = '\t')[,3]

plot(covData, type= 'l', xlab = 'Position', ylab = 'Depth', col = 'blue')
title(paste('Coverage depth for sample', nameinfo[[1]][3]), outer=TRUE)

dev.off()
