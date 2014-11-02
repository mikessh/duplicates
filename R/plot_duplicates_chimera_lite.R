require(ggplot2)

args<-commandArgs(TRUE)

prefix<-args[1]

df <- read.table(paste(prefix, "ilmn_dupl_stats.txt", sep = "_"), header = TRUE, sep = "\t")

pdf(paste(prefix, "ilmn_dupl_stats.pdf", sep = "_"))
ggplot() + geom_density(data = subset(df, sample == "same_sample"), aes(dist), alpha=.3, adjust = 5) + 
scale_y_continuous("Density") + 
scale_x_log10("Cluster distance, pixels") +
geom_vline(xintercept=c(100), linetype="dotted") +
facet_grid(umi ~ sample + cdr3)
dev.off()