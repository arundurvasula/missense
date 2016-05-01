library(data.table)
library(MASS)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

par(mfrow=c(2,3))

## MAF > 0.1
system("awk 'FNR==1 && NR!=1 { while (/^gen/) getline; } 1 {print}' results/selection/maf-0.1*.param > results/selection/pooled-maf-0.1.param")
a <- fread("results/selection/pooled-maf-0.1.param")

## MAF < 0.1
system("awk 'FNR==1 && NR!=1 { while (/^gen/) getline; } 1 {print}' results/selection/maf-lt-0.1*.param > results/selection/pooled-maf-lt-0.1.param")
b <- fread("results/selection/pooled-maf-lt-0.1.param")

## Syn
system("awk 'FNR==1 && NR!=1 { while (/^gen/) getline; } 1 {print}' results/selection/syn*.param > results/selection/pooled-syn.param")
c <- fread("results/selection/pooled-syn.param")

## Plot it
pdf("./selection-age.pdf",width=8.5,height=11)
image(kde2d(b$alpha1/(2*3500), b$alpha2/(2*3500), n=200), col=r, main="Selection coefficients for\nrare missense mutations", xlab=expression("s"[2]), ylab=expression("s"[1]))
image(kde2d(a$alpha1/(2*3500), a$alpha2/(2*3500), n=200), col=r, main="Selection coefficients for\ncommon missense mutations", xlab=expression("s"[2]), ylab=expression("s"[1]))
image(kde2d(c$alpha1/(2*3500), c$alpha2/(2*3500), n=200), col=r, main="Selection coefficients for\nsynonymous mutations", xlab=expression("s"[2]), ylab=expression("s"[1]))
plot(density(b$age*2*3500), col="orange", main="Age distribution of\nrare missense mutations", xlab="Allele age (years)")
plot(density(a$age*2*3500), col="orange", main="Age distribution of\ncommon missense mutations", xlab="Allele age (years)")
plot(density(c$age*2*3500), col="orange", main="Age distribution of\nsynonymous mutations", xlab="Allele age (years)")
dev.off()
