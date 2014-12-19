# Read in data from allele_model.R
mcmc_samples <- readRDS("model2.Rdata")

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