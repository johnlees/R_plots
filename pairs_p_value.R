# 429 mutations observed total
# 2000kb genome => 429/2000 mutations per kb
# 2000 genes to correct for multiple testing

df <- as.data.frame(matrix(0,ncol=3,nrow=11*3))
i <- 0
for (gene_len in c(1,4,10)) {
  for (mutations in seq(0,10)) {
    df[i*11 + mutations+1,3] = paste(gene_len,"kb",sep='')
    df[i*11 + mutations+1,2] <- as.numeric(p.adjust(1-sum(dpois(seq(0,mutations),429/2000*gene_len)), method="holm", n=2000))
    df[i*11 + mutations+1,1] = mutations
  }
  i <- i+1
}
colnames(df) <- c("mutations", "p_value","operon_length")

ggplot(df, aes(x=mutations,y=p_value,colour=factor(operon_length))) + geom_point() + theme_bw(base_size=18)  + scale_x_discrete(limits=seq(0,10,1),breaks=seq(0,10,1)) + scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"),name="operon length",breaks=c("1kb","4kb","10kb"),labels=c("1kb","4kb","10kb")) + geom_line() + geom_hline(yintercept=0.05,colour="red")