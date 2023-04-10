############################################
### effects of corrosion inhibitor levels
### on microbial community diversity
### Lou LaMartina, finalized Apr 10, 2023
############################################

# test for normality, print results
shapiro.verbose <- function(x) {
  stat = shapiro.test(x)$p.value
  cat("p =", stat)
  if (stat >= 0.05) { 
    cat(" - not normal distribution - use wilcox.test\n")  
  } else { cat(" - normal distribution - use t.test\n") }
}


setwd("~/Desktop/Postdoc/Kimbell/GitHub")
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(reshape2)
library(ape)
library(ggrepel)



#################
### load data ###
#################
# removed A2_Day0_2 - low quality

# ASV counts; change hyphen to underscore
counts <- read.csv("Data/Microcosm_ASV_counts.csv", stringsAsFactors = F)
counts$Sample_name <- gsub("-", "_", counts$Sample_name)
rownames(counts) <- counts$Sample_name
counts <- counts[-1]


# sample metadata
info <- read.csv("Data/Microcosm_sample_info.csv", stringsAsFactors = F)
rownames(info) <- info$Sample_ID
info <- info[order(info$File_name),]


# change row names to sample ID
identical(info$File_name, rownames(counts))
rownames(counts) <- rownames(info)


# convert to relative abundance
relabun <- counts / rowSums(counts)




######################
### beta diversity ###
######################

########
### stat

# stat
beta <- data.frame(ID1 = rownames(relabun), scores(vegdist(relabun, "bray")))


# add info
beta <- melt(beta, variable.name = "ID2", value.name = "Bray", id.vars = "ID1")
beta$ID2 <- as.character(beta$ID2)
beta$Compare <- NA
for (i in 1:nrow(beta)) {
  beta$Compare[i] <- paste(sort(unlist(beta[i, c("ID1", "ID2")])), 
                             collapse = "__")
}
beta <- unique(beta[beta$Bray > 0, c("Compare", "Bray")])
beta$ID1 <- sapply(strsplit(beta$Compare, "__"), '[', 1)
beta$ID2 <- sapply(strsplit(beta$Compare, "__"), '[', 2)
beta$Set1 <- sapply(strsplit(beta$ID1, "_"), '[', 1)
beta$Set2 <- sapply(strsplit(beta$ID2, "_"), '[', 1)
beta$Treat1 <- sapply(strsplit(beta$ID1, "_"), '[', 2)
beta$Treat2 <- sapply(strsplit(beta$ID2, "_"), '[', 2)
beta$Day1 <- sapply(strsplit(beta$ID1, "_"), '[', 3)
beta$Day2 <- sapply(strsplit(beta$ID2, "_"), '[', 3)



############
### organize

# untreated water, normal (set2) vs high (set3) inhibitor concentrations
# compared same day, same set, to same day & set controls
beta <- rbind(
  beta[beta$Set1 == "set3" & beta$Set2 == "set3" &
           beta$Treat1 == "control" & 
           beta$Day1 == beta$Day2, -1],
  beta[beta$Set1 == "set2" & beta$Set2 == "set2" &
           beta$Treat1 == "control" & 
           beta$Day1 == beta$Day2, -1])


# don't want day zero; comparing treatments to controls on same day
beta <- beta[beta$Day1 != "day0",]
beta$Treat3 <- beta$Treat2
beta$Treat3[beta$Treat3 != "control"] <- "treatment"


# summary stats for plot
beta.stats <- data.frame(aggregate(Bray ~ Set1 + Treat1 + Treat3, mean, data = beta),
                               aggregate(Bray ~ Set1 + Treat1 + Treat3, sd, data = beta)[4])
colnames(beta.stats)[4:5] <- c("mean", "sd")



###################
### normality tests

# normal levels - beta-curtis btwn replicates of same-day controls
normctrl.beta <- beta[beta$Set1 == "set2" & beta$Treat2 == "control",]
cat(round(mean(normctrl.beta$Bray), digits = 2), "±", round(sd(normctrl.beta$Bray), digits = 2))
shapiro.verbose(normctrl.beta$Bray) 
# 0.21 ± 0.04
# p = 0.5197797 - not normal distribution - use wilcox.test


# normal levels - beta-curtis btwn replicate inhibitors and same-day controls
normexp.beta <- beta[beta$Set1 == "set2" & beta$Treat2 != "control",]
cat(round(mean(normexp.beta$Bray), digits = 2), "±", round(sd(normexp.beta$Bray), digits = 2))
shapiro.verbose(normexp.beta$Bray) 
# 0.41 ± 0.16
# p = 1.674899e-05 - not normal distribution - use wilcox.test


# high levels - beta-curtis btwn replicates of same-day controls
hictrl.beta <- beta[beta$Set1 == "set3" & beta$Treat2 == "control",]
cat(round(mean(hictrl.beta$Bray), digits = 2), "±", round(sd(hictrl.beta$Bray), digits = 2))
shapiro.verbose(hictrl.beta$Bray) 
# 0.27 ± 0.06
# p = 0.9751937 - not normal distribution - use wilcox.test


# high levels - beta-curtis btwn replicate inhibitors and same-day controls
hiexp.beta <- beta[beta$Set1 == "set3" & beta$Treat2 != "control",]
cat(round(mean(hiexp.beta$Bray), digits = 2), "±", round(sd(hiexp.beta$Bray), digits = 2))
shapiro.verbose(hiexp.beta$Bray) 
# 0.54 ± 0.2
# p = 3.580506e-06 - normal distribution - use t.test



###############
### comparisons

# normal inhibitor levels - 
# were exp samples a greater distance from same-day controls
# then controls to replicate controls?
cat("p =", wilcox.test(x = normexp.beta$Bray, y = normctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 0.0002936866


# high inhibitor levels - 
# were exp samples a greater distance from same-day controls
# then controls to replicate controls?
cat("p =", wilcox.test(x = hiexp.beta$Bray, y = hictrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 0.002784969



########
### plot

beta.plot <-
  ggplot(beta, aes(x = Treat3, y = Bray)) +
  
  # points
  geom_jitter(width = 0.1, shape = 1, size = 0.5) +
  geom_errorbar(data = beta.stats,
                aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, linewidth = 1) +
  geom_point(data = beta.stats, aes(y = mean), size = 3) +
  
  # significance
  geom_text(data = data.frame(Treat3 = 1.5, Bray = 0.83, Set1 = "set2"), 
            label = "p = 0.00029", size = 3) +
  geom_text(data = data.frame(Treat3 = 1.5, Bray = 0.83, Set1 = "set3"), 
            label = "p = 0.0028", size = 3) +
  annotate("segment", x = 1, xend = 2, y = 0.8, yend = 0.8) +
  annotate("segment", x = c(1,2), xend = c(1,2), y = 0.8, yend = 0.78) +
  
  # labels
  facet_grid(~ Set1, labeller = labeller(Set1 = c(set2 = "Normal inhibitor levels", 
                                                  set3 = "High inhibitor levels"))) +
  scale_x_discrete(labels = c(control = "No inhibitor\n(control)", 
                              treatment = "With inhibitor\n(experiment)")) + 
  
  # theme
  theme_classic() +
  theme(axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 10, color = "black", face = "bold"),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(linewidth = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25)) +
  guides(color = guide_legend(keyheight = 1.5)) +
  labs(x = "Treatment", y = "Bray-Curtis dissimilarity score\ncompared to same-day control")
beta.plot
ggsave("../Plots/beta_inhibitorLevels.pdf", plot = beta.plot, device = "pdf", width = 6, height = 4, units = "in")




#######################
### alpha diversity ###
#######################

########
### stat 

# stat
alpha <- data.frame(Sample_ID = rownames(relabun), Shannon = diversity(relabun, "shannon"))


# add info
alpha <- merge(alpha, info, by = "Sample_ID")
alpha <- alpha[alpha$Set != "set1",]


# summary stats for plot
alpha.stats <- data.frame(aggregate(Shannon ~ Set + Treatment, mean, data = alpha),
                          aggregate(Shannon ~ Set + Treatment, sd, data = alpha)[3])
colnames(alpha.stats)[3:4] <- c("mean", "sd")



###################
### normality tests

# normal levels - shannon diversity in controls
normctrl.alpha <- alpha[alpha$Set == "set2" & alpha$Treatment == "control",]
cat(round(mean(normctrl.alpha$Shannon), digits = 2), "±", round(sd(normctrl.alpha$Shannon), digits = 2))
shapiro.verbose(normctrl.alpha$Shannon) 
# 5.37 ± 0.36
# p = 0.2625465 - not normal distribution - use wilcox.test


# normal levels - shannon diversity in experiments
normexp.alpha <- alpha[alpha$Set == "set2" & alpha$Treatment != "control",]
cat(round(mean(normexp.alpha$Shannon), digits = 2), "±", round(sd(normexp.alpha$Shannon), digits = 2))
shapiro.verbose(normexp.alpha$Shannon) 
# 5.46 ± 0.31
# p = 0.07917687 - not normal distribution - use wilcox.test


# hi levels - shannon diversity in controls
hictrl.alpha <- alpha[alpha$Set == "set3" & alpha$Treatment == "control",]
cat(round(mean(hictrl.alpha$Shannon), digits = 2), "±", round(sd(hictrl.alpha$Shannon), digits = 2))
shapiro.verbose(hictrl.alpha$Shannon) 
# 5.07 ± 0.32
# p = 0.1758912 - not normal distribution - use wilcox.test


# hi levels - shannon diversity in experiments
hiexp.alpha <- alpha[alpha$Set == "set3" & alpha$Treatment != "control",]
cat(round(mean(hiexp.alpha$Shannon), digits = 2), "±", round(sd(hiexp.alpha$Shannon), digits = 2))
shapiro.verbose(hiexp.alpha$Shannon) 
# 4.9 ± 0.34
# p = 0.2871963 - not normal distribution - use wilcox.test



###############
### comparisons

# normal inhibitor levels - 
# differences in community richness between amended samples and controls?
cat("p =", wilcox.test(x = normctrl.alpha$Shannon, y = normexp.alpha$Shannon, paired = F)$p.value)
# p = 0.4613567


# high inhibitor levels - 
# were exp samples a greater distance from same-day controls
# then controls to replicate controls?
cat("p =", wilcox.test(x = hiexp.alpha$Shannon, y = hictrl.alpha$Shannon, paired = F)$p.value)
# p = 0.1450496




########
### plot

alpha.plot <-
  ggplot(alpha, aes(x = Treatment, y = Shannon)) +
  
  # points
  geom_jitter(width = 0.1, shape = 1, size = 0.5) +
  geom_errorbar(data = alpha.stats,
                aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, linewidth = 1) +
  geom_point(data = alpha.stats, aes(y = mean), size = 3) +
  
  # significance
  geom_text(data = data.frame(Treatment = 1.5, Shannon = 6.1, Set = "set2"), 
            label = "p = 0.46", size = 3) +
  geom_text(data = data.frame(Treatment = 1.5, Shannon = 6.1, Set = "set3"), 
            label = "p = 0.15", size = 3) +
  annotate("segment", x = 1, xend = 2, y = 6, yend = 6) +
  annotate("segment", x = c(1,2), xend = c(1,2), y = 6, yend = 5.95) +
  
  # labels
  facet_grid(~ Set, labeller = labeller(Set = c(set2 = "Normal inhibitor levels", 
                                                set3 = "High inhibitor levels"))) +
  scale_x_discrete(labels = c(control = "No inhibitor\n(control)", 
                              treatment = "With inhibitor\n(experiment)")) + 
  
  # theme
  theme_classic() +
  theme(axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 10, color = "black", face = "bold"),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(linewidth = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25)) +
  guides(color = guide_legend(keyheight = 1.5)) +
  labs(x = "Treatment", y = "Shannon alpha diversity score")
alpha.plot
ggsave("../Plots/alpha_inhibitorLevels.pdf", plot = alpha.plot, device = "pdf", width = 6, height = 4, units = "in")
