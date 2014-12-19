library(ggplot2)
library(plyr)

#Here, V2 is X, V3 is Y and V1 is binned giving three boxes at each X
#(use a factor e.g.
 odds_bin <- rep("log(OR) < 2",length(power_meningitis$V1))
 odds_bin[abs(log(power_meningitis$V1)) > 2] <- "2 < log(OR) < 4"
 odds_bin[abs(log(power_meningitis$V1)) > 4] <- "log(OR) > 4"

ggplot(power_meningitis) + geom_boxplot(aes(x=factor(round_any(V2,0.1)),y=V3,fill=factor(odds_bin,levels = levels(factor(odds_bin))[c(3,1,2)]))) + xlab('MAF') + ylab('Samples') + labs(fill='log(OR)') + theme_set(theme_gray(base_size=14)) + labs(title="80% power for 35% cases")

#White background, 18pt font
#theme_bw(base_size=18)