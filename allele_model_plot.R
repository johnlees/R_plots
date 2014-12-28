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
pairs_file <- paste(data_location, "strep_pairs.txt", sep="/")

# Outcomes
alleles = c("A", "B", "C", "D", "E", "F")
tissues = c("CSF", "Blood")

# Top level parameters
parameters = c("mu[1,1]", "mu[2,1]", "mu[1,2]", "mu[2,2]", "mu[1,3]", "mu[2,3]", 
	"mu[1,4]", "mu[2,4]", "mu[1,5]", "mu[2,5]", "mu[1,6]", "mu[2,6]")

#
# Functions
#
intervals <- function(chains,file_out="confidence_intervals.txt")
{
	conf_intervals <- p.interval(chains,HPD=TRUE)
	means <- colMeans(chains)

	capture.output(conf_intervals,file_out)	
	return(means, conf_intervals)
}


#
# Main
#

#
# Read in data
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
gelman <- gelman.diag(mcmc_samples)
heidel <- heidel.diag(mcmc_samples)
capture.output(list("gelman"=gelman,"heidel"=heidel),file="convergence_diagnostics.txt")

pdf("kappa_gelman.pdf")
gelman.plot(mcmc_samples[,"kappa",drop=FALSE])
dev.off()

# Plot of kappa + chain, ggplot of kappa
pdf("kappa_chain.pdf")
plot(mcmc_samples[,"kappa",drop=F])
dev.off()

pdf("kappa_posterior.pdf")
ggplot(kappa_data, aes(x=kappa)) 
+ geom_histogram(aes(y=..density..), binwidth=0.025,colour="black",fill="white") 
+ geom_density(alpha=0.2, fill="#FF9999")
dev.off()

# means and 95% HPD intervals, concatenating all chains together
(conf_intervals, means) = intervals(all_chains,"confidence_intervals.txt")

# 6*2 means and 95% CIs for mus
# 	geom pointrange plot of this
#   could also try plotting multiple geom_density, for each allele
tissue_alleles <- as.data.frame(rep(alleles,length(tissues)))
tissue_alleles$tissue <- rep(tissues,length(alleles))
tissue_alleles$means = means[parameters]
tissue_alleles$upper = conf_intervals[parameters,"Upper"]
tissue_alleles$lower = conf_intervals[parameters,"Lower"]

ggplot(tissue_alleles) + 
geom_pointrange(aes(x=alleles,y=means,ymax=upper,ymin=lower,colour=tissue),position=position_dodge(width=0.9)) 
+ xlab("Allele") + ylab("Proportion in tissue") + ylim(0,1)

# 6*1192 means and 95% CIs for pis
# 6*(1192/2) means and 95% CIs for pi[csf] - pi[blood] for all pairs
for (i in 1:nrow(paired_samples))
{
	csf_sample = paste("pi[",paired_samples$CSF_index[i],",",sep='')
	blood_sample = paste("pi[",paired_samples$Blood_index[i],",",sep='')

	for (j in 1:length(alleles))
	{
		csf_allele = paste(csf_sample,j,"]",sep='')
		blood_allele = paste(blood_sample,j,"]",sep='')
		
		mcmc_diff[,i+j] = all_chains[,csf_sample] - all_chains[,blood_sample]
		row_headers = c(row_headers, paste(paired_samples$samples[i],"_",alleles[j])
	}
}

rownames(mcmc_diff) = row_headers

(pairs_conf_intervals, pairs_means) = intervals(mcmc_diff,"pairs_intervals.txt")
