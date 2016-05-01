library(data.table)
library(MASS)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

par(mfrow=c(2,3))

## MAF > 0.1
system("awk 'FNR==1 && NR!=1 { while (/^gen/) getline; } 1 {print}' results/selection/maf-0.1*.param > results/selection/pooled-maf-0.1.param")
a <- fread("results/selection/pooled-maf-0.1.param")
image(kde2d(a$alpha1/(2*3500), a$alpha2/(2*3500), n=200), col=r, main="Selection coefficients for common missense mutations", xlab=expression("s"[2]), ylab=expression("s"[1]))
plot(density(a$age*2*3500), col="orange", main="Age distribution of common missense mutations", xlab="Allele age (years)")

## MAF < 0.1
system("cat results/selection/maf-lt-0.1*.param > results/selection/pooled-maf-lt-0.1.param")
b <- fread("results/selection/pooled-maf-lt-0.1.param")
image(kde2d(b$alpha1/(2*3500), b$alpha2/(2*3500), n=200), col=r, main="Selection coefficients for rare missense mutations", xlab=expression("s"[2]), ylab=expression("s"[1]))
plot(density(b$age*2*3500), col="orange", main="Age distribution of rare missense mutations", xlab="Allele age (years)")

## Syn
system("cat results/selection/syn*.param > results/selection/pooled-syn.param")
c <- fread("results/selection/pooled-syn.param")
image(kde2d(c$alpha1/(2*3500), c$alpha2/(2*3500), n=200), col=r, main="Selection coefficients for synonymous mutations", xlab=expression("s"[2]), ylab=expression("s"[1]))
plot(density(c$age*2*3500), col="orange", main="Age distribution of synonymous mutations", xlab="Allele age (years)")


