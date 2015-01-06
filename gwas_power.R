library(ggplot2)
library(plyr)

power_meningitis <- read.delim("~/Documents/PhD/GWAS/power_simulation/fixed_sample_size_power.txt", header=F)
colnames(power_meningitis) = c("OR","MAF","Power")

#Here, V2 is X, V3 is Y and V1 is binned giving three boxes at each X
#(use a factor e.g.
 odds_bin <- rep("log(OR) < 2",length(power_meningitis$OR))
 odds_bin[abs(log(power_meningitis$OR)) > 2] <- "2 < log(OR) < 4"
 odds_bin[abs(log(power_meningitis$OR)) > 4] <- "log(OR) > 4"

# Fixed power, variable sample size
ggplot(power_meningitis) + 
  geom_boxplot(aes(x=factor(round_any(V2,0.1)),y=V3,fill=factor(odds_bin,levels = levels(factor(odds_bin))[c(3,1,2)]))) + 
  xlab('MAF') + ylab('Samples') + labs(fill='log(OR)') + theme_set(theme_gray(base_size=14)) + 
  labs(title="80% power for 35% cases")

# Fixed sample size, variable power
ggplot(power_meningitis) + 
  geom_boxplot(aes(x=factor(round_any(MAF,0.1)),y=Power,fill=factor(odds_bin,levels = levels(factor(odds_bin))[c(3,1,2)]))) + 
  xlab('MAF') + ylab('Power') + labs(fill='log(OR)') + theme_set(theme_gray(base_size=14)) + 
  labs(title="Power for 35% cases, 1200 samples")

#White background, 18pt font
#theme_bw(base_size=18)