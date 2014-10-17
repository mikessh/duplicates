require(ggplot2)

args<-commandArgs(TRUE)

prefix<-args[0]

df <- read.table(paste(prefix, "ilmn_dupl_stats.txt", sep = "_"), header = TRUE)

# Manually melt to subsets
df.d <- subset(df, dist < 100)
df.q <- data.frame(dist = df$dist, qual = df$qual, set = "Same")
df.q <- rbind(df.q, data.frame(dist = df.d$dist, qual = df.d$qual, set = "Same-close"))
df.q <- rbind(df.q, data.frame(dist = df$bg.dist, qual = df$bg.qual, set = "Distinct"))

# Will plot using facets
df.qq <- df.q
df.q$facet <- rep("Distribution", nrow(df.q))
df.qq$facet <- rep("By quality", nrow(df.q))
df.q <- rbind(df.q, df.qq)

pdf(paste(prefix, "ilmn_dupl_stats.pdf", sep = "_"))

ggplot()+
  geom_density(data=subset(df.q, facet=="Distribution" & set %in% c("Same", "Distinct")), aes(x = dist, fill = set), alpha=.3)+
  stat_density2d(data=subset(df.q,facet=="By quality"), aes(x = dist, y = qual, color = set, n = 200)) +
  geom_vline(xintercept=c(100), linetype="dotted") +
  scale_y_continuous("") +
  scale_x_log10("Cluster distance, pixels") +
  facet_grid(facet~.,scales="free_y") +
  guides(fill = guide_legend("Read set", order = 1), color = guide_legend("", order = 2))

dev.off()