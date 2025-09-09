# First get column of f(BBAA)/total

infile<-read.table("Calonectris.nohybrids_Calonectris.filtered.maxmiss75.mac3.nohybrids_Dmin.pop_pairs_analysis.csv", header=T)
Ber<-infile[infile$P1=="Ber" | infile$P2 =="Ber",]
Ber <- Ber[order(Ber$f.BBAA.),]
write.table(Ber, "Ber_sorted_fBBAA.csv", row.names=F, quote=F)
