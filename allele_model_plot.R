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

# Data input
data_location <- "~/Documents/PhD/hsd_locus/mapping"

mcmc_file <- paste(data_location, "model2.Rdata", sep="/")
pairs_file <- paste(data_location, "pair_rows.txt", sep="/")

# Outcomes
alleles = c("A", "B", "C", "D", "E", "F")
tissues = c("CSF", "Blood")

# Top level parameters
parameters = c("mu[1,1]", "mu[2,1]", "mu[1,2]", "mu[2,2]", "mu[1,3]", "mu[2,3]", 
	"mu[1,4]", "mu[2,4]", "mu[1,5]", "mu[2,5]", "mu[1,6]", "mu[2,6]")

#
# Main
#

#
# mcmc.list data structure
#
# Sample means - pi
# First index is sample number
# Second index is allele: 1 = A; 2 = B ... 6 = E
# summary(coda_samples2[,"pi[1,1]",drop=FALSE])
# plot(coda_samples2[,"pi[1,1]",drop=FALSE])
#
# Tissue means - mu
# First index is tissue: 1 = csf; 2 = blood
# Second index is allele: 1 = A; 2 = B ... 6 = E
# summary(coda_samples2[,"mu[1,2]",drop=FALSE])
# plot(coda_samples2[,"mu[1,2]",drop=FALSE])
#
# Relation between mu and alpha - kappa
# summary(coda_samples2[,"kappa",drop=FALSE])
# plot(coda_samples2[,"kappa",drop=FALSE])
#
# Priors for sample level Dirichlet - alpha
# Indicies as for my
# summary(coda_samples2[,"alpha[1,2]",drop=FALSE])
# plot(coda_samples2[,"alpha[1,2]",drop=FALSE])

#
# Read in data
#

# Read in data from allele_model.R
mcmc_samples <- readRDS(mcmc_file)
all_chains <- as.matrix(mcmc_samples)

# Pair data
paired_samples <- read.delim(pairs_file)

# Chain diagnostics plots
# Autocorrelation
pdf("kappa_autocorr.pdf")
autocorr.plot(mcmc_samples[,c("kappa","mu[1,1]","mu[2,1]"),drop=FALSE])
dev.off()

# Chain convergence
gelman <- gelman.diag(mcmc_samples[,c(parameters,"kappa"),drop=FALSE])
heidel <- heidel.diag(mcmc_samples[,c(parameters,"kappa"),drop=FALSE])
capture.output(list("gelman"=gelman,"heidel"=heidel),file="convergence_diagnostics.txt")

pdf("kappa_gelman.pdf")
gelman.plot(mcmc_samples[,"kappa",drop=FALSE])
dev.off()

# Plot of kappa + chain, ggplot of kappa
pdf("kappa_chain.pdf")
plot(mcmc_samples[,"kappa",drop=F])
dev.off()

kappa_data = as.data.frame(all_chains[,"kappa"])
colnames(kappa_data) = "kappa"

pdf("kappa_posterior.pdf")
ggplot(kappa_data, aes(x=kappa)) + 
geom_histogram(aes(y=..density..), binwidth=0.025,colour="black",fill="white") + 
geom_density(alpha=0.2, fill="#FF9999")
dev.off()

# means and 95% HPD intervals, concatenating all chains together
conf_intervals <- p.interval(all_chains,HPD=TRUE)
capture.output(conf_intervals,file="confidence_intervals.txt")
means = colMeans(all_chains)

# 6*2 means and 95% CIs for mus
# 	geom pointrange plot of this
#   could also try plotting multiple geom_density, for each allele
tissue_alleles <- as.data.frame(rep(alleles,each=length(tissues)))
colnames(tissue_alleles)[1] = "alleles"
tissue_alleles$tissue <- rep(tissues,length(alleles))
tissue_alleles$means = means[parameters]
tissue_alleles$upper = conf_intervals[parameters,"Upper"]
tissue_alleles$lower = conf_intervals[parameters,"Lower"]

pdf("tissue_alleles.pdf",width=10,height=7)
ggplot(tissue_alleles) + 
geom_pointrange(aes(x=alleles,y=means,ymax=upper,ymin=lower,colour=tissue),position=position_dodge(width=0.9),size=0.8) + 
xlab("Allele") + ylab("Proportion in tissue") + ylim(0,1) + theme_grey(base_size = 18)
dev.off()

# 6*1192 means and 95% CIs for pis
# 6*(1192/2) means and 95% CIs for pi[csf] - pi[blood] for all pairs
mcmc_diff <- matrix(data=NA,nrow=nrow(all_chains),ncol=(nrow(paired_samples)*length(alleles)))
col_headers <- vector(mode="character",length=(nrow(paired_samples)*length(alleles)))

for (i in 1:nrow(paired_samples))
{
	csf_sample = paste("pi[",paired_samples$CSF_index[i],",",sep='')
	blood_sample = paste("pi[",paired_samples$Blood_index[i],",",sep='')

	for (j in 1:length(alleles))
	{
		column_number = (i-1)*length(alleles) + j
    
    csf_allele = paste(csf_sample,j,"]",sep='')
		blood_allele = paste(blood_sample,j,"]",sep='')
		
		mcmc_diff[,column_number] = all_chains[,csf_allele] - all_chains[,blood_allele]
    col_headers[column_number] = paste(paired_samples$Sample[i],alleles[j],sep='_')
	}
}
colnames(mcmc_diff) = col_headers

# Print out confidence intervals
pair_conf_intervals <- p.interval(mcmc_diff,HPD=TRUE)
capture.output(pair_conf_intervals,file="pairs_intervals.txt")