#
# Libraries
#
library(coda)
library(LaplacesDemon)
library(ggplot2)
library(plyr)

#
# Constants
#

# Outcomes
alleles = c("A", "B", "C", "D", "E", "F")
tissues = c("CSF", "Blood")

# Top level parameters
parameters = c("mu[1,1]", "mu[2,1]", "mu[1,2]", "mu[2,2]", "mu[1,3]", "mu[2,3]", 
	"mu[1,4]", "mu[2,4]", "mu[1,5]", "mu[2,5]", "mu[1,6]", "mu[2,6]")

# Read in data from allele_model.R
mcmc_samples <- readRDS("model2.Rdata")
all_chains <- as.matrix(mcmc_samples)

# Chain diagnostics plots
# Autocorrelation
pdf("kappa_autocorr.pdf")
autocorr.plot(mcmc_samples[,c("kappa","mu[1,1]",drop=FALSE])
dev.off()

# Chain convergence
gelman <- gelman.diag(mcmc_samples)
heidel <- heidel.diag(mcmc_samples)

pdf("kappa_gelman.pdf")
gelman.plot(mcmc_samples[,"kappa",drop=FALSE])
dev.off()

# Plot of kappa + chain, ggplot of kappa
pdf("kappa_chain.pdf")
plot(mcmc_samples[,"kappa",drop=F])
dev.off()

pdf("kappa_density.pdf")
kappa_data <- as.data.frame(all_chains[,"kappa"])
names(kappa_data) <- "kappa"
# TODO: add dashes to bottom of this
ggplot(kappa_data) + geom_density(aes(x=kappa))
dev.off()
# TODO: plot as histogram too
# ggplot(kappa_data) + geom_histogram(aes(x=kappa,fill='red'))

# means and 95% HPD intervals, concatenating all chains together
conf_intervals <- p.interval(all_chains,HPD=TRUE)
means <- colMeans(all_chains)

# 6*2 means and 95% CIs for mus
# 	plot of this
#   or should it just be a boxplot?
tissue_alleles <- as.data.frame(c(alleles,alleles))
tissue_alleles$means = means[c("mu[1]","mu[2]")]
# etc

names(tissue_alleles) <- c("allele", "tissue", "mean", "upper", "lower")
for (i in 1:length(alleles) by 2)
{
	for (j in 1:length(tissues))
	{
		index = i + j - 1
		tissue_alleles[index,] = c(alleles[index], tissues[j], means[parameters[index]], conf_intervals[parameters[index],])	
	}
}

ggplot(tissue_alleles) + geom_pointrange(aes(x=c("N terminus"),y=means,ymax=upper,ymin=lower,colour=tissue),position=position_dodge(width=0.9)) + xlab("N terminus") + ylab("theta") + ylim(0,1)


# 6*1192 means and 95% CIs for pis
# 6*(1192/2) means and 95% CIs for pi[csf] - pi[blood] for all pairs

# Will also need to read in pairs

# Sample means - pi
# First index is sample number
# Second index is allele: 1 = A; 2 = B ... 6 = E
summary(coda_samples2[,"pi[1,1]",drop=FALSE])
plot(coda_samples2[,"pi[1,1]",drop=FALSE])

# Tissue means - mu
# First index is tissue: 1 = csf; 2 = blood
# Second index is allele: 1 = A; 2 = B ... 6 = E
summary(coda_samples2[,"mu[1,2]",drop=FALSE])
plot(coda_samples2[,"mu[1,2]",drop=FALSE])

# Relation between mu and alpha - kappa
summary(coda_samples2[,"kappa",drop=FALSE])
plot(coda_samples2[,"kappa",drop=FALSE])

# Priors for sample level Dirichlet - alpha
# Indicies as for my
summary(coda_samples2[,"alpha[1,2]",drop=FALSE])
plot(coda_samples2[,"alpha[1,2]",drop=FALSE])

#
# Extracting mean and 95% conf
# (remember to do as.matrix only once to save a lot of time) 
#
mean(as.matrix(coda_samples1[,c("kappa"),drop=FALSE]))
quantile(as.matrix(coda_samples1[,c("kappa"),drop=FALSE]),0.025)
quantile(as.matrix(coda_samples1[,c("kappa"),drop=FALSE]),0.975)

# better
library(LaplacesDemon)
p.interval(as.matrix(coda_samples[,"kappa",drop=F]),HPD=TRUE)