#!/bin/Rscript

library(ggplot2)

infile <- read.table("Calonectris.CALBOR.cv.error", header=F, col.names=c("K","CVE"))
a<-ggplot(infile, aes(x=K, y=CVE)) + geom_line(alpha=0.5, , linewidth= 1.2, colour="darkblue") + geom_point(size=3.5, colour="darkblue") + 
    theme_light() + scale_x_continuous (breaks=seq(0,10,1))
ggsave("Calonectris.CALBOR.cv.error.pdf", a, width=4, height=4)

infile <- read.table("Calonectris.CALDIO.cv.error", header=F, col.names=c("K","CVE"))
a<-ggplot(infile, aes(x=K, y=CVE)) + geom_line(alpha=0.5, , linewidth= 1.2, colour="darkorange") + geom_point(size=3.5, colour="darkorange") + theme_light() + scale_x_continuous (breaks=seq(0,10,1))
ggsave("Calonectris.CALDIO.cv.error.pdf", a, width=4, height=4)

infile <- read.table("Calonectris.noEDW.cv.error", header=F, col.names=c("K","CVE"))
a<-ggplot(infile, aes(x=K, y=CVE)) + geom_line(alpha=0.5, , linewidth= 1.2, colour="firebrick3") + geom_point(size=3.5, colour="firebrick3") + theme_light() + scale_x_continuous (breaks=seq(0,10,1))
ggsave("Calonectris.noEDW.cv.error.pdf", a, width=4, height=4)

infile <- read.table("Calonectris.EDW.cv.error", header=F, col.names=c("K","CVE"))
a<-ggplot(infile, aes(x=K, y=CVE)) + geom_line(alpha=0.5, , linewidth= 1.2, colour="forestgreen") + geom_point(size=3.5, colour="forestgreen") + theme_light() + scale_x_continuous (breaks=seq(0,10,1))
ggsave("Calonectris.EDW.cv.error.pdf", a, width=4, height=4)
ggsave("Calonectris.EDW.cv.error.jpg", a, width=4, height=4)

infile <- read.table("Calonectris.LEU.cv.error", header=F, col.names=c("K","CVE"))
infile<-infile[infile$K>0.5,]      # Error with K10 when running ADMIXTURE
a<-ggplot(infile, aes(x=K, y=CVE)) + geom_line(alpha=0.5, , linewidth= 1.2, colour="grey50") + geom_point(size=3.5, colour="grey50") + theme_light() + scale_x_continuous (breaks=seq(0,10,1))
ggsave("Calonectris.LEU.cv.error.pdf", a, width=4, height=4)
