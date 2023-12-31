BYneg=read.table("ICL1-GFP_20231006_BY_negctrl2_028.txt") 

A13.4h=read.table("ICL1-GFP_20231006_A1-3_4h2_023.txt")
A13.24h=read.table("ICL1-GFP_20231006_A1-3_24h2_029.txt")
A110.4h=read.table("ICL1-GFP_20231006_A1-10_4h2_015.txt")
A110.24h=read.table("ICL1-GFP_20231006_A1-10_24h2_030.txt")
A23.4h=read.table("ICL1-GFP_20231006_A2-3_4h2_016.txt")
A23.24h=read.table("ICL1-GFP_20231006_A2-3_24h2_031.txt")
A25.4h=read.table("ICL1-GFP_20231006_A2-5_4h2_017.txt")
A25.24h=read.table("ICL1-GFP_20231006_A2-5_24h2_032.txt")
A31.4h=read.table("ICL1-GFP_20231006_A3-1_4h2_018.txt")
A31.24h=read.table("ICL1-GFP_20231006_A3-1_24h2_033.txt")
A36.4h=read.table("ICL1-GFP_20231006_A3-6_4h2_019.txt")
A36.24h=read.table("ICL1-GFP_20231006_A3-6_24h2_034.txt")

A21.4h=read.table("PDC1-GFP_20231009_A2-1_ICL1-4h_028.txt")
A21.24h=read.table("PDC1-GFP_20231009_A2-1_ICL1-24h_027.txt")

N11.4h=read.table("ICL1-GFP_20231006_N1-1_4h2_022.txt")
N11.24h=read.table("ICL1-GFP_20231006_N1-1_24h2_035.txt")
N112.4h=read.table("ICL1-GFP_20231006_N1-12_4h2_021.txt")
N112.24h=read.table("ICL1-GFP_20231006_N1-12_24h2_036.txt")
N22.4h=read.table("ICL1-GFP_20231006_N2-2_4h2_024.txt")
N22.24h=read.table("ICL1-GFP_20231006_N2-2_24h2_037.txt")
N26.4h=read.table("ICL1-GFP_20231006_N2-6_4h2_025.txt")
N26.24h=read.table("ICL1-GFP_20231006_N2-6_24h2_038.txt")
N37.4h=read.table("ICL1-GFP_20231006_N3-7_4h2_026.txt")
N37.24h=read.table("ICL1-GFP_20231006_N3-7_24h2_039.txt")
N39.4h=read.table("ICL1-GFP_20231006_N3-9_4h2_027.txt")
N39.24h=read.table("ICL1-GFP_20231006_N3-9_24h2_040.txt")


## Example FSC cutoff 

plot(A13.4h$FSC.A/A13.4h$FSC.H,A13.4h$FSC.A,pch=20,col=rgb(0,0,0,0.05),xlim=c(0.5,5),log="y", xlab="FSC-A/FSC-H",ylab="FSC-A",cex.axis=1.2,cex.lab=1.4)
abline(v=1.9,lwd=2,lty=2,col="red")
dev.copy2pdf(file="A13.4h_example_cutoff_FSC.pdf",useDingbats=F)
dev.off()


#plot(A13.4h$SSC.A/A13.4h$SSC.H,A13.4h$SSC.A,pch=20,col=rgb(0,0,0,0.05),xlim=c(0.5,5),log="y")
# plot(A13.4h$FSC.A/A13.4h$FSC.H,A13.4h$SSC.A/A13.4h$SSC.H,pch=20,col=rgb(0,0,0,0.05),xlim=c(0.5,5),log="y")


## Plot the distribution of FITC vs FSC

x1=log(A13.4h$FSC.A[which(A13.4h$FSC.A>2000 & A13.4h$SSC.A>500 & A13.4h$FITC.A>0 & A13.4h$FSC.A/A13.4h$FSC.H >= 1 & A13.4h$FSC.A/A13.4h$FSC.H <= 1.9)])
x2=log(A13.4h$FITC.A[which(A13.4h$FSC.A>2000 & A13.4h$SSC.A>500 & A13.4h$FITC.A>0 & A13.4h$FSC.A/A13.4h$FSC.H >= 1 & A13.4h$FSC.A/A13.4h$FSC.H <= 1.9)])


df<-data.frame(x1,x2)
x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),heights=c(1,1))
par(mar = c(4.5,4.5,2.5,1), font = 2)
plot(x2~x1, data=df[order(df$dens),], pch=16, col=col, cex=1, xlim=c(6.5,11.5), ylim=c(0,12),xlab="Cell size (a.u.)",ylab="GFP expression level (a.u.)",cex.lab=1.4,cex.axis=1.5)

## abline(v=7.3,lwd=2,lty=2)
par(mar = c(3.5,3,3,2.5))

min=1
max=0
for(i in 1:256)
{
  x=which(df$dens==i)
  d=length(x)/length(x1)
  if(d<min) { min=d }
  if(d>max) { max=d }
}
if(min==0) { min=1/length(x1) }

ColorLevels=seq(min,max,length=length(cols))
image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),col=cols,xlab="",ylab="",xaxt="n", las = 1,cex.axis=1.3)

# Save the image file in pdf format

dev.copy2pdf(file="A13.4h_FSC_vs_FITC.pdf",useDingbats=F)
dev.off()

rm(df)
rm(x)
rm(x1)
rm(x2)


x1=log(N11.4h$FSC.A[which(N11.4h$FSC.A>2000 & N11.4h$SSC.A>500 & N11.4h$FITC.A>0 & N11.4h$FSC.A/N11.4h$FSC.H >= 1 & N11.4h$FSC.A/N11.4h$FSC.H <= 1.9)])
x2=log(N11.4h$FITC.A[which(N11.4h$FSC.A>2000 & N11.4h$SSC.A>500 & N11.4h$FITC.A>0 & N11.4h$FSC.A/N11.4h$FSC.H >= 1 & N11.4h$FSC.A/N11.4h$FSC.H <= 1.9)])


df<-data.frame(x1,x2)
x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),heights=c(1,1))
par(mar = c(4.5,4.5,2.5,1), font = 2)
plot(x2~x1, data=df[order(df$dens),], pch=16, col=col, cex=1, xlim=c(6.5,11.5), ylim=c(0,12),xlab="Cell size (a.u.)",ylab="GFP expression level (a.u.)",cex.lab=1.4,cex.axis=1.5)

## abline(v=7.3,lwd=2,lty=2)
par(mar = c(3.5,3,3,2.5))

min=1
max=0
for(i in 1:256)
{
  x=which(df$dens==i)
  d=length(x)/length(x1)
  if(d<min) { min=d }
  if(d>max) { max=d }
}
if(min==0) { min=1/length(x1) }

ColorLevels=seq(min,max,length=length(cols))
image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),col=cols,xlab="",ylab="",xaxt="n", las = 1,cex.axis=1.3)

# Save the image file in pdf format

dev.copy2pdf(file="N11.4h_FSC_vs_FITC.pdf",useDingbats=F)
dev.off()

rm(df)
rm(x)
rm(x1)
rm(x2)


## Boxplot / Vioplot 

by.final=BYneg$FITC.A[which(BYneg$FSC.A>2000 & BYneg$SSC.A>500 & BYneg$FITC.A>0 & BYneg$FSC.A/BYneg$FSC.H >= 1 & BYneg$FSC.A/BYneg$FSC.H <= 1.9)]

A13.4h.final=A13.4h$FITC.A[which(A13.4h$FSC.A>2000 & A13.4h$SSC.A>500 & A13.4h$FITC.A>0 & A13.4h$FSC.A/A13.4h$FSC.H >= 1 & A13.4h$FSC.A/A13.4h$FSC.H <= 1.9)]
A13.24h.final=A13.24h$FITC.A[which(A13.24h$FSC.A>2000 & A13.24h$SSC.A>500 & A13.24h$FITC.A>0 & A13.24h$FSC.A/A13.24h$FSC.H >= 1 & A13.24h$FSC.A/A13.24h$FSC.H <= 1.9)]

A110.4h.final=A110.4h$FITC.A[which(A110.4h$FSC.A>2000 & A110.4h$SSC.A>500 & A110.4h$FITC.A>0 & A110.4h$FSC.A/A110.4h$FSC.H >= 1 & A110.4h$FSC.A/A110.4h$FSC.H <= 1.9)]
A110.24h.final=A110.24h$FITC.A[which(A110.24h$FSC.A>2000 & A110.24h$SSC.A>500 & A110.24h$FITC.A>0 & A110.24h$FSC.A/A110.24h$FSC.H >= 1 & A110.24h$FSC.A/A110.24h$FSC.H <= 1.9)]

#A23.4h.final=A23.4h$FITC.A[which(A23.4h$FSC.A>2000 & A23.4h$SSC.A>500 & A23.4h$FITC.A>0 & A23.4h$FSC.A/A23.4h$FSC.H >= 1 & A23.4h$FSC.A/A23.4h$FSC.H <= 1.9)]
#A23.24h.final=A23.24h$FITC.A[which(A23.24h$FSC.A>2000 & A23.24h$SSC.A>500 & A23.24h$FITC.A>0 & A23.24h$FSC.A/A23.24h$FSC.H >= 1 & A23.24h$FSC.A/A23.24h$FSC.H <= 1.9)]

A21.4h.final=A21.4h$FITC.A[which(A21.4h$FSC.A>2000 & A21.4h$SSC.A>500 & A21.4h$FITC.A>0 & A21.4h$FSC.A/A21.4h$FSC.H >= 1 & A21.4h$FSC.A/A21.4h$FSC.H <= 1.9)]
A21.24h.final=A21.24h$FITC.A[which(A21.24h$FSC.A>2000 & A21.24h$SSC.A>500 & A21.24h$FITC.A>0 & A21.24h$FSC.A/A21.24h$FSC.H >= 1 & A21.24h$FSC.A/A21.24h$FSC.H <= 1.9)]

A25.4h.final=A25.4h$FITC.A[which(A25.4h$FSC.A>2000 & A25.4h$SSC.A>500 & A25.4h$FITC.A>0 & A25.4h$FSC.A/A25.4h$FSC.H >= 1 & A25.4h$FSC.A/A25.4h$FSC.H <= 1.9)]
A25.24h.final=A25.24h$FITC.A[which(A25.24h$FSC.A>2000 & A25.24h$SSC.A>500 & A25.24h$FITC.A>0 & A25.24h$FSC.A/A25.24h$FSC.H >= 1 & A25.24h$FSC.A/A25.24h$FSC.H <= 1.9)]

A31.4h.final=A31.4h$FITC.A[which(A31.4h$FSC.A>2000 & A31.4h$SSC.A>500 & A31.4h$FITC.A>0 & A31.4h$FSC.A/A31.4h$FSC.H >= 1 & A31.4h$FSC.A/A31.4h$FSC.H <= 1.9)]
A31.24h.final=A31.24h$FITC.A[which(A31.24h$FSC.A>2000 & A31.24h$SSC.A>500 & A31.24h$FITC.A>0 & A31.24h$FSC.A/A31.24h$FSC.H >= 1 & A31.24h$FSC.A/A31.24h$FSC.H <= 1.9)]

A36.4h.final=A36.4h$FITC.A[which(A36.4h$FSC.A>2000 & A36.4h$SSC.A>500 & A36.4h$FITC.A>0 & A36.4h$FSC.A/A36.4h$FSC.H >= 1 & A36.4h$FSC.A/A36.4h$FSC.H <= 1.9)]
A36.24h.final=A36.24h$FITC.A[which(A36.24h$FSC.A>2000 & A36.24h$SSC.A>500 & A36.24h$FITC.A>0 & A36.24h$FSC.A/A36.24h$FSC.H >= 1 & A36.24h$FSC.A/A36.24h$FSC.H <= 1.9)]

N11.4h.final=N11.4h$FITC.A[which(N11.4h$FSC.A>2000 & N11.4h$SSC.A>500 & N11.4h$FITC.A>0 & N11.4h$FSC.A/N11.4h$FSC.H >= 1 & N11.4h$FSC.A/N11.4h$FSC.H <= 1.9)]
N11.24h.final=N11.24h$FITC.A[which(N11.24h$FSC.A>2000 & N11.24h$SSC.A>500 & N11.24h$FITC.A>0 & N11.24h$FSC.A/N11.24h$FSC.H >= 1 & N11.24h$FSC.A/N11.24h$FSC.H <= 1.9)]

N112.4h.final=N112.4h$FITC.A[which(N112.4h$FSC.A>2000 & N112.4h$SSC.A>500 & N112.4h$FITC.A>0 & N112.4h$FSC.A/N112.4h$FSC.H >= 1 & N112.4h$FSC.A/N112.4h$FSC.H <= 1.9)]
N112.24h.final=N112.24h$FITC.A[which(N112.24h$FSC.A>2000 & N112.24h$SSC.A>500 & N112.24h$FITC.A>0 & N112.24h$FSC.A/N112.24h$FSC.H >= 1 & N112.24h$FSC.A/N112.24h$FSC.H <= 1.9)]

N22.4h.final=N22.4h$FITC.A[which(N22.4h$FSC.A>2000 & N22.4h$SSC.A>500 & N22.4h$FITC.A>0 & N22.4h$FSC.A/N22.4h$FSC.H >= 1 & N22.4h$FSC.A/N22.4h$FSC.H <= 1.9)]
N22.24h.final=N22.24h$FITC.A[which(N22.24h$FSC.A>2000 & N22.24h$SSC.A>500 & N22.24h$FITC.A>0 & N22.24h$FSC.A/N22.24h$FSC.H >= 1 & N22.24h$FSC.A/N22.24h$FSC.H <= 1.9)]

N26.4h.final=N26.4h$FITC.A[which(N26.4h$FSC.A>2000 & N26.4h$SSC.A>500 & N26.4h$FITC.A>0 & N26.4h$FSC.A/N26.4h$FSC.H >= 1 & N26.4h$FSC.A/N26.4h$FSC.H <= 1.9)]
N26.24h.final=N26.24h$FITC.A[which(N26.24h$FSC.A>2000 & N26.24h$SSC.A>500 & N26.24h$FITC.A>0 & N26.24h$FSC.A/N26.24h$FSC.H >= 1 & N26.24h$FSC.A/N26.24h$FSC.H <= 1.9)]

N37.4h.final=N37.4h$FITC.A[which(N37.4h$FSC.A>2000 & N37.4h$SSC.A>500 & N37.4h$FITC.A>0 & N37.4h$FSC.A/N37.4h$FSC.H >= 1 & N37.4h$FSC.A/N37.4h$FSC.H <= 1.9)]
N37.24h.final=N37.24h$FITC.A[which(N37.24h$FSC.A>2000 & N37.24h$SSC.A>500 & N37.24h$FITC.A>0 & N37.24h$FSC.A/N37.24h$FSC.H >= 1 & N37.24h$FSC.A/N37.24h$FSC.H <= 1.9)]

N39.4h.final=N39.4h$FITC.A[which(N39.4h$FSC.A>2000 & N39.4h$SSC.A>500 & N39.4h$FITC.A>0 & N39.4h$FSC.A/N39.4h$FSC.H >= 1 & N39.4h$FSC.A/N39.4h$FSC.H <= 1.9)]
N39.24h.final=N39.24h$FITC.A[which(N39.24h$FSC.A>2000 & N39.24h$SSC.A>500 & N39.24h$FITC.A>0 & N39.24h$FSC.A/N39.24h$FSC.H >= 1 & N39.24h$FSC.A/N39.24h$FSC.H <= 1.9)]


##boxplot(by.final,A13.4h.final,A110.4h.final,A21.4h.final,A25.4h.final,A31.4h.final,A36.4h.final,N11.4h.final,N112.4h.final,N22.4h.final,N26.4h.final,N37.4h.final,N39.4h.final,pch=16,ylab="ICL1 expression (a.u.)",cex.lab=1.4,log="y")
##axis(side=1, at=c(1,2,3,4,5,6,7), labels=c("NC","A1","A2","A3","N1","N2","N3"), las=2, cex.axis=0.6)
##dev.copy2pdf(file="Boxplot_A&N_ICL1GFP_4h_20231006.pdf",useDingbats=F)

=========================

## Density plot 4h 

plot(density(A13.4h.final),xlim=c(0,1000),ylim=c(0,0.008),lwd=2,col="#000066",cex.axis=1.2,cex.lab=1.4,xlab="ICL1 expression (a.u.)", main="ICL1 expression")
lines(density(A110.4h.final),lwd=2,col="#000066")
lines(density(A21.4h.final),lwd=2,col="#3333FF")
lines(density(A25.4h.final),lwd=2,col="#3333FF")
lines(density(A31.4h.final),lwd=2,col="#3399FF")
lines(density(A36.4h.final),lwd=2,col="#3399FF")

lines(density(N11.4h.final),lwd=2,col="#CC3300")
lines(density(N112.4h.final),lwd=2,col="#CC3300")
lines(density(N22.4h.final),lwd=2,col="#FF3300")
lines(density(N26.4h.final),lwd=2,col="#FF3300")
lines(density(N37.4h.final),lwd=2,col="#FF9900")
lines(density(N39.4h.final),lwd=2,col="#FF9900")

lines(density(by.final),lwd=2,col="grey")


abline(v=mlv(A13.4h.final,method="venter"),lty=2,col="#000066")
abline(v=mlv(A110.4h.final,method="venter"),lty=2,col="#000066")
abline(v=mlv(A21.4h.final,method="venter"),lty=2,col="#3333FF")
abline(v=mlv(A25.4h.final,method="venter"),lty=2,col="#3333FF")
abline(v=mlv(A31.4h.final,method="venter"),lty=2,col="#3399FF")
abline(v=mlv(A36.4h.final,method="venter"),lty=2,col="#3399FF")
abline(v=mlv(N11.4h.final,method="venter"),lty=2,col="#CC3300")
abline(v=mlv(N112.4h.final,method="venter"),lty=2,col="#CC3300")
abline(v=mlv(N22.4h.final,method="venter"),lty=2,col="#FF3300")
abline(v=mlv(N26.4h.final,method="venter"),lty=2,col="#FF3300")
abline(v=mlv(N37.4h.final,method="venter"),lty=2,col="#FF9900")
abline(v=mlv(N39.4h.final,method="venter"),lty=2,col="#FF9900")
abline(v=mlv(by.final,method="venter"),lty=2,col="darkgrey")

dev.copy2pdf(file="ICL1_4h_density.pdf", width=9, height=6, useDingbats=F)
dev.off()


## Density plot 24h 

plot(density(A13.24h.final),xlim=c(0,12000),ylim=c(0,0.0005),lwd=2,col="#000066",cex.axis=1.2,cex.lab=1.4,xlab="ICL1 expression (a.u.)", main="ICL1 expression")
lines(density(A110.24h.final),lwd=2,col="#000066")
lines(density(A21.24h.final),lwd=2,col="#3333FF")
lines(density(A25.24h.final),lwd=2,col="#3333FF")
lines(density(A31.24h.final),lwd=2,col="#3399FF")
lines(density(A36.24h.final),lwd=2,col="#3399FF")

lines(density(N11.24h.final),lwd=2,col="#CC3300")
lines(density(N112.24h.final),lwd=2,col="#CC3300")
lines(density(N22.24h.final),lwd=2,col="#FF3300")
lines(density(N26.24h.final),lwd=2,col="#FF3300")
lines(density(N37.24h.final),lwd=2,col="#FF9900")
lines(density(N39.24h.final),lwd=2,col="#FF9900")

lines(density(by.final),lwd=2,col="grey")


abline(v=mlv(A13.24h.final,method="venter"),lty=2,col="#000066")
abline(v=mlv(A110.24h.final,method="venter"),lty=2,col="#000066")
abline(v=mlv(A21.24h.final,method="venter"),lty=2,col="#3333FF")
abline(v=mlv(A25.24h.final,method="venter"),lty=2,col="#3333FF")
abline(v=mlv(A31.24h.final,method="venter"),lty=2,col="#3399FF")
abline(v=mlv(A36.24h.final,method="venter"),lty=2,col="#3399FF")
abline(v=mlv(N11.24h.final,method="venter"),lty=2,col="#CC3300")
abline(v=mlv(N112.24h.final,method="venter"),lty=2,col="#CC3300")
abline(v=mlv(N22.24h.final,method="venter"),lty=2,col="#FF3300")
abline(v=mlv(N26.24h.final,method="venter"),lty=2,col="#FF3300")
abline(v=mlv(N37.24h.final,method="venter"),lty=2,col="#FF9900")
abline(v=mlv(N39.24h.final,method="venter"),lty=2,col="#FF9900")
abline(v=mlv(by.final,method="venter"),lty=2,col="darkgrey")

dev.copy2pdf(file="ICL1_24h_density.pdf", width=9, height=6, useDingbats=F)
dev.off()


===================

## Histogram 4h

par(mfrow=c(3,4),omi=c(0.1,0.25,0.1,0.05),mar=c(1,1,1,1))

par(mar=c(6,2,1.5,2))

hist(prob=T,A13.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#000066",border="#000066",xlab="",ylab="",main="")
abline(v=mlv(A13.4h.final,method="venter"),lwd=2,lty=2,col="red")

hist(prob=T,A110.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#000066",border="#000066",xlab="",ylab="",main="")
abline(v=mlv(A110.4h.final,method="venter"),lwd=2, lty=2,col="red")

hist(prob=T,N11.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#CC3300",border="#CC3300",xlab="",ylab="",main="")
abline(v=mlv(N11.4h.final,method="venter"),lwd=2,lty=2,col="blue")

hist(prob=T,N112.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#CC3300",border="#CC3300",xlab="",ylab="",main="")
abline(v=mlv(N112.4h.final,method="venter"),lwd=2,lty=2,col="blue")


hist(prob=T,A21.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#3333FF",border="#3333FF",xlab="",ylab="",main="")
abline(v=mlv(A21.4h.final,method="venter"),lwd=2,lty=2,col="red")

hist(prob=T,A25.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#3333FF",border="#3333FF",xlab="",ylab="",main="")
abline(v=mlv(A25.4h.final,method="venter"),lwd=2,lty=2,col="red")


hist(prob=T,N22.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#FF3300",border="#FF3300",xlab="",ylab="",main="")
abline(v=mlv(N22.4h.final,method="venter"),lwd=2,lty=2,col="blue")

hist(prob=T,N26.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#FF3300",border="#FF3300",xlab="",ylab="",main="")
abline(v=mlv(N26.4h.final,method="venter"),lwd=2,lty=2,col="blue")


hist(prob=T,A31.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#3399FF",border="#3399FF",xlab="",ylab="",main="")
abline(v=mlv(A31.4h.final,method="venter"),lwd=2,lty=2,col="red")

hist(prob=T,A36.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#3399FF",border="#3399FF",xlab="",ylab="",main="")
abline(v=mlv(A36.4h.final,method="venter"),lwd=2,lty=2,col="red")


hist(prob=T,N37.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#FF9900",border="#FF9900",xlab="",ylab="",main="")
abline(v=mlv(N37.4h.final,method="venter"),lwd=2,lty=2,col="blue")

hist(prob=T,N39.4h.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="#FF9900",border="#FF9900",xlab="",ylab="",main="")
abline(v=mlv(N39.4h.final,method="venter"),lwd=2,lty=2,col="blue")


hist(prob=T,by.final,breaks=2000,xlim=c(0,1000),ylim=c(0,0.004),col="grey",border="grey",xlab="",ylab="",main="")
abline(v=mlv(by.final,method="venter"),lwd=2,lty=2,col="darkgrey")


dev.copy2pdf(file="ICL1_4h_histogram.pdf", width=9, height=6, useDingbats=F)
dev.off()



## Histogram 24h

par(mfrow=c(3,4),omi=c(0.1,0.25,0.1,0.05),mar=c(1,1,1,1))

par(mar=c(6,2,1.5,2))

hist(prob=T,A13.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#000066",border="#000066",xlab="",ylab="",main="")
#abline(v=mlv(A13.24h.final,method="venter"),lwd=2,lty=2,col="red")

hist(prob=T,A110.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#000066",border="#000066",xlab="",ylab="",main="")
#abline(v=mlv(A110.24h.final,method="venter"),lwd=2, lty=2,col="red")

hist(prob=T,N11.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#CC3300",border="#CC3300",xlab="",ylab="",main="")
#abline(v=mlv(N11.24h.final,method="venter"),lwd=2,lty=2,col="blue")

hist(prob=T,N112.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#CC3300",border="#CC3300",xlab="",ylab="",main="")
#abline(v=mlv(N112.24h.final,method="venter"),lwd=2,lty=2,col="blue")


hist(prob=T,A21.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#3333FF",border="#3333FF",xlab="",ylab="",main="")
#abline(v=mlv(A21.24h.final,method="venter"),lwd=2,lty=2,col="red")

hist(prob=T,A25.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#3333FF",border="#3333FF",xlab="",ylab="",main="")
#abline(v=mlv(A25.24h.final,method="venter"),lwd=2,lty=2,col="red")

hist(prob=T,N22.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#FF3300",border="#FF3300",xlab="",ylab="",main="")
#abline(v=mlv(N22.24h.final,method="venter"),lwd=2,lty=2,col="blue")

hist(prob=T,N26.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#FF3300",border="#FF3300",xlab="",ylab="",main="")
#abline(v=mlv(N26.24h.final,method="venter"),lwd=2,lty=2,col="blue")

hist(prob=T,A31.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#3399FF",border="#3399FF",xlab="",ylab="",main="")
#abline(v=mlv(A31.24h.final,method="venter"),lwd=2,lty=2,col="red")

hist(prob=T,A36.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#3399FF",border="#3399FF",xlab="",ylab="",main="")
#abline(v=mlv(A36.24h.final,method="venter"),lwd=2,lty=2,col="red")

hist(prob=T,N37.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#FF9900",border="#FF9900",xlab="",ylab="",main="")
#abline(v=mlv(N37.24h.final,method="venter"),lwd=2,lty=2,col="blue")

hist(prob=T,N39.24h.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="#FF9900",border="#FF9900",xlab="",ylab="",main="")
#abline(v=mlv(N39.24h.final,method="venter"),lwd=2,lty=2,col="blue")

hist(prob=T,by.final,breaks=2000,xlim=c(0,12000),ylim=c(0,0.0004),col="grey",border="grey",xlab="",ylab="",main="")
#abline(v=mlv(by.final,method="venter"),lwd=2,lty=2,col="darkgrey")


dev.copy2pdf(file="ICL1_24h_histogram.pdf", width=9, height=6, useDingbats=F)
dev.off()
