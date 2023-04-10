#########################################
### effects of corrosion inhibitor types
### on microbial community diversity
### Lou LaMartina, finalized Apr 10, 2023
#########################################

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
beta.stats <- data.frame(aggregate(Bray ~ Set1 + Treat2, mean, data = beta),
                         aggregate(Bray ~ Set1 + Treat2, sd, data = beta)[3])
colnames(beta.stats)[3:4] <- c("mean", "sd")




###############
### overall ###
###############

###################
### normality tests

# beta-curtis btwn replicates of same-day controls
ctrl.beta <- beta[beta$Treat2 == "control",]
cat(round(mean(ctrl.beta$Bray), digits = 2), "±", round(sd(ctrl.beta$Bray), digits = 2))
shapiro.verbose(ctrl.beta$Bray)
# 0.24 ± 0.06
# p = 0.9381506 - not normal distribution - use wilcox.test


# beta-curtis btwn Na2SiO3 samples and same-day controls
ss.beta <- beta[beta$Treat2 == "Na2SiO3",]
cat(round(mean(ss.beta$Bray), digits = 2), "±", round(sd(ss.beta$Bray), digits = 2))
shapiro.verbose(ss.beta$Bray)
# 0.47 ± 0.22
# p = 1.312329e-05 - normal distribution - use t.test


# beta-curtis btwn NaPO4 samples and same-day controls
so.beta <- beta[beta$Treat2 == "NaPO4",]
cat(round(mean(so.beta$Bray), digits = 2), "±", round(sd(so.beta$Bray), digits = 2))
shapiro.verbose(so.beta$Bray)
# 0.31 ± 0.07
# p = 0.3138244 - not normal distribution - use wilcox.test


# beta-curtis btwn NaPO4 samples and same-day controls
zo.beta <- beta[beta$Treat2 == "ZnPO4",]
cat(round(mean(zo.beta$Bray), digits = 2), "±", round(sd(zo.beta$Bray), digits = 2))
shapiro.verbose(zo.beta$Bray)
# 0.65 ± 0.07
# p = 0.002127465 - normal distribution - use t.test



###############
### comparisons

# were distances btwn Na2SiO3-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = ss.beta$Bray, y = ctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 0.0002856006


# were distances btwn NaPO4-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = so.beta$Bray, y = ctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 0.0003515537


# were distances btwn ZnPO4-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = zo.beta$Bray, y = ctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 1.435368e-11




##############
### normal ###
##############

###################
### normality tests

# beta-curtis btwn replicates of same-day controls
normctrl.beta <- beta[beta$Treat2 == "control" & beta$Set2 == "set2",]
cat(round(mean(normctrl.beta$Bray), digits = 2), "±", round(sd(normctrl.beta$Bray), digits = 2))
shapiro.verbose(normctrl.beta$Bray) 
# 0.21 ± 0.04
# p = 0.5197797 - not normal distribution - use wilcox.test


# beta-curtis btwn Na2SiO3 samples and same-day controls
normss.beta <- beta[beta$Treat2 == "Na2SiO3" & beta$Set2 == "set2",]
cat(round(mean(normss.beta$Bray), digits = 2), "±", round(sd(normss.beta$Bray), digits = 2))
shapiro.verbose(normss.beta$Bray) 
# 0.26 ± 0.04
# p = 0.8423096 - not normal distribution - use wilcox.test


# beta-curtis btwn NaPO4 samples and same-day controls
normso.beta <- beta[beta$Treat2 == "NaPO4" & beta$Set2 == "set2",]
cat(round(mean(normso.beta$Bray), digits = 2), "±", round(sd(normso.beta$Bray), digits = 2))
shapiro.verbose(normso.beta$Bray) 
# 0.36 ± 0.06
# p = 0.822262 - not normal distribution - use wilcox.test


# beta-curtis btwn NaPO4 samples and same-day controls
normzo.beta <- beta[beta$Treat2 == "ZnPO4" & beta$Set2 == "set2",]
cat(round(mean(normzo.beta$Bray), digits = 2), "±", round(sd(normzo.beta$Bray), digits = 2))
shapiro.verbose(normzo.beta$Bray) 
# 0.63 ± 0.02
# p = 0.7176429 - not normal distribution - use wilcox.test



###############
### comparisons

# were distances btwn normal levels of Na2SiO3-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = normss.beta$Bray, y = normctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 0.01121876


# were distances btwn normal levels of NaPO4-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = normso.beta$Bray, y = normctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 2.971857e-05


# were distances btwn normal levels of ZnPO4-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = normzo.beta$Bray, y = normctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 7.429641e-06




############
### high ###
############

###################
### normality tests

# beta-curtis btwn replicates of same-day controls
hictrl.beta <- beta[beta$Treat2 == "control" & beta$Set2 == "set3",]
cat(round(mean(hictrl.beta$Bray), digits = 2), "±", round(sd(hictrl.beta$Bray), digits = 2))
shapiro.verbose(hictrl.beta$Bray) 
# 0.27 ± 0.06
# p = 0.9751937 - not normal distribution - use wilcox.test


# beta-curtis btwn Na2SiO3 samples and same-day controls
hiss.beta <- beta[beta$Treat2 == "Na2SiO3" & beta$Set2 == "set3",]
cat(round(mean(hiss.beta$Bray), digits = 2), "±", round(sd(hiss.beta$Bray), digits = 2))
shapiro.verbose(hiss.beta$Bray) 
# 0.69 ± 0.05
# p = 0.2578805 - not normal distribution - use wilcox.test


# beta-curtis btwn NaPO4 samples and same-day controls
hiso.beta <- beta[beta$Treat2 == "NaPO4" & beta$Set2 == "set3",]
cat(round(mean(hiso.beta$Bray), digits = 2), "±", round(sd(hiso.beta$Bray), digits = 2))
shapiro.verbose(hiso.beta$Bray) 
# 0.27 ± 0.04
# p = 0.3456629 - not normal distribution - use wilcox.test


# beta-curtis btwn NaPO4 samples and same-day controls
hizo.beta <- beta[beta$Treat2 == "ZnPO4" & beta$Set2 == "set3",]
cat(round(mean(hizo.beta$Bray), digits = 2), "±", round(sd(hizo.beta$Bray), digits = 2))
shapiro.verbose(hizo.beta$Bray) 
# 0.67 ± 0.09
# p = 0.005740826 - normal distribution - use t.test



###############
### comparisons

# were distances btwn hial levels of Na2SiO3-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = hiss.beta$Bray, y = hictrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# 7.429641e-06


# were distances btwn hial levels of NaPO4-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = hiso.beta$Bray, y = hictrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 0.3851303


# were distances btwn hial levels of ZnPO4-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = hizo.beta$Bray, y = hictrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 7.429641e-06


# NaPO4 more of an impact at normal levels than high?
cat("p =", wilcox.test(x = normso.beta$Bray, y = hiso.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 1.719236e-05


# differences btwn ZnPO4?
cat("p =", wilcox.test(x = normzo.beta$Bray, y = hizo.beta$Bray, paired = F)$p.value)
# p = 0.5010106


# Na2SiO3 more of an impact at normal levels than high?
cat("p =", wilcox.test(x = normss.beta$Bray, y = hiss.beta$Bray, paired = F, alternative = "less")$p.value)
# p = 1.719236e-05


# differences btwn controls?
cat("p =", wilcox.test(x = normctrl.beta$Bray, y = hictrl.beta$Bray, paired = F)$p.value)
# p = 0.09307359



########
### plot

inhiblabs <- c(Na2SiO3 = expression(Na[2]*SiO[3]),
               NaPO4 = expression(NaPO[4]),
               ZnPO4 = expression(ZnPO[4]))


beta.plot <-
  ggplot(beta, aes(x = Treat2, y = Bray, color = Set1)) +
  
  # points
  geom_jitter(shape = 1, size = 0.5, position = 
                position_jitterdodge(jitter.width = 0.1, dodge.width = 0.25)) +
  geom_errorbar(data = beta.stats,
                aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, linewidth = 1, position = position_dodge(0.25)) +
  geom_point(data = beta.stats, aes(y = mean), size = 3, position = position_dodge(0.25)) +
  
  scale_x_discrete(labels = inhiblabs) +
  scale_color_manual(labels = c(set2 = "normal", set3 = "high"),
                       values = brewer.pal(8, "Set1")[2:1]) +
  
  # theme
  theme_classic() +
  theme(axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        legend.position = c(0.15,0.85),
        legend.background = element_rect(fill = NA, color = "grey80"),
        legend.title = element_text(size = 9, color = "black", face = "bold"),
        strip.text = element_text(size = 10, color = "black", face = "bold"),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(linewidth = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25)) +
  guides(color = guide_legend(keyheight = 1)) +
  labs(x = "Corrosion inhibitor type", color = "Inhibitor\nlevel",
       y = "Bray-Curtis dissimilarity score\ncompared to same-day control")
beta.plot
ggsave("../Plots/beta_inhibitorTypes.pdf", plot = beta.plot, device = "pdf", width = 4, height = 4, units = "in")




#######################
### alpha diversity ###
#######################

########
### stat 

# stat
alpha <- data.frame(Sample_ID = rownames(relabun), Shannon = diversity(relabun, "shannon"))


# add info
alpha <- merge(alpha, info, by = "Sample_ID")
alpha <- alpha[alpha$Set != "set1" & alpha$Day != "day0",]


# summary stats for plot
alpha.stats <- data.frame(aggregate(Shannon ~ Set + Inhibitor, mean, data = alpha),
                          aggregate(Shannon ~ Set + Inhibitor, sd, data = alpha)[3])
colnames(alpha.stats)[3:4] <- c("mean", "sd")



###################
### normality tests

### normal levels

# shannon diversity in controls
normctrl.alpha <- alpha[alpha$Set == "set2" & alpha$Treatment == "control",]
cat(round(mean(normctrl.alpha$Shannon), digits = 2), "±", round(sd(normctrl.alpha$Shannon), digits = 2))
shapiro.verbose(normctrl.alpha$Shannon) 
# 5.54 ± 0.21
# p = 0.305246 - not normal distribution - use wilcox.test


# shannon diversity in Na2SiO3
normss.alpha <- alpha[alpha$Set == "set2" & alpha$Inhibitor == "Na2SiO3",]
cat(round(mean(normss.alpha$Shannon), digits = 2), "±", round(sd(normss.alpha$Shannon), digits = 2))
shapiro.verbose(normss.alpha$Shannon) 
# 5.52 ± 0.34
# p = 0.3759599 - not normal distribution - use wilcox.test


# shannon diversity in NaPO4
normso.alpha <- alpha[alpha$Set == "set2" & alpha$Inhibitor == "NaPO4",]
cat(round(mean(normso.alpha$Shannon), digits = 2), "±", round(sd(normso.alpha$Shannon), digits = 2))
shapiro.verbose(normso.alpha$Shannon) 
# 5.66 ± 0.18
# p = 0.07938339 - not normal distribution - use wilcox.test


# shannon diversity in ZnPO4
normzo.alpha <- alpha[alpha$Set == "set2" & alpha$Inhibitor == "ZnPO4",]
cat(round(mean(normzo.alpha$Shannon), digits = 2), "±", round(sd(normzo.alpha$Shannon), digits = 2))
shapiro.verbose(normzo.alpha$Shannon) 
# 5.2 ± 0.23
# p = 0.08383906 - not normal distribution - use wilcox.test


### high levels 

# shannon diversity in controls
hictrl.alpha <- alpha[alpha$Set == "set3" & alpha$Treatment == "control",]
cat(round(mean(hictrl.alpha$Shannon), digits = 2), "±", round(sd(hictrl.alpha$Shannon), digits = 2))
shapiro.verbose(hictrl.alpha$Shannon) 
# 5.15 ± 0.35
# p = 0.1313504 - not normal distribution - use wilcox.test


# high levels - shannon diversity in Na2SiO3
hiss.alpha <- alpha[alpha$Set == "set3" & alpha$Inhibitor == "Na2SiO3",]
cat(round(mean(hiss.alpha$Shannon), digits = 2), "±", round(sd(hiss.alpha$Shannon), digits = 2))
shapiro.verbose(hiss.alpha$Shannon) 
# 4.73 ± 0.25
# p = 0.38999 - not normal distribution - use wilcox.test


# high levels - shannon diversity in NaPO4
hiso.alpha <- alpha[alpha$Set == "set3" & alpha$Inhibitor == "NaPO4",]
cat(round(mean(hiso.alpha$Shannon), digits = 2), "±", round(sd(hiso.alpha$Shannon), digits = 2))
shapiro.verbose(hiso.alpha$Shannon) 
# 5.25 ± 0.29
# p = 0.3540485 - not hial distribution - use wilcox.test


# high levels - shannon diversity in ZnPO4
hizo.alpha <- alpha[alpha$Set == "set3" & alpha$Inhibitor == "ZnPO4",]
cat(round(mean(hizo.alpha$Shannon), digits = 2), "±", round(sd(hizo.alpha$Shannon), digits = 2))
shapiro.verbose(hizo.alpha$Shannon) 
# 4.71 ± 0.14
# p = 0.6359253 - not normal distribution - use wilcox.test



###############
### comparisons

### normal inhibitor levels

# differences in community richness between Na2SiO3 samples and controls?
cat("p =", wilcox.test(x = normctrl.alpha$Shannon, y = normss.alpha$Shannon, paired = F)$p.value)
# p = 0.9372294


# differences in community richness between NaPO4 samples and controls?
cat("p =", wilcox.test(x = normctrl.alpha$Shannon, y = normso.alpha$Shannon, paired = F)$p.value)
# p = 0.2402597


# normal inhibitor levels - 
# differences in community richness between ZnPO4 samples and controls?
cat("p =", wilcox.test(x = normctrl.alpha$Shannon, y = normzo.alpha$Shannon, paired = F, alternative = "greater")$p.value)
# p = 0.02056277


### high inhibitor levels

# differences in community richness between Na2SiO3 samples and controls?
cat("p =", wilcox.test(x = hictrl.alpha$Shannon, y = hiss.alpha$Shannon, paired = F, alternative = "greater")$p.value)
# p = 0.02056277


# differences in community richness between NaPO4 samples and controls?
cat("p =", wilcox.test(x = hictrl.alpha$Shannon, y = hiso.alpha$Shannon, paired = F)$p.value)
# p = 0.8181818


# differences in community richness between ZnPO4 samples and controls?
cat("p =", wilcox.test(x = hictrl.alpha$Shannon, y = hizo.alpha$Shannon, paired = F, alternative = "greater")$p.value)
# p = 0.02056277



########
### plot

alpha.plot <-
  ggplot(alpha, aes(x = Inhibitor, y = Shannon, color = Set)) +
  
  # points
  geom_jitter(shape = 1, size = 0.5, position = 
                position_jitterdodge(jitter.width = 0.1, dodge.width = 0.25)) +
  geom_errorbar(data = alpha.stats,
                aes(y = mean, ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, linewidth = 1, position = position_dodge(0.25)) +
  geom_point(data = alpha.stats, aes(y = mean), size = 3, position = position_dodge(0.25)) +
  
  scale_x_discrete(labels = inhiblabs) +
  scale_color_manual(labels = c(set2 = "normal", set3 = "high"),
                     values = brewer.pal(8, "Set1")[2:1]) +
  
  # theme
  theme_classic() +
  theme(axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        legend.position = c(0.15,0.85),
        legend.background = element_rect(fill = NA, color = "grey80"),
        legend.title = element_text(size = 9, color = "black", face = "bold"),
        strip.text = element_text(size = 10, color = "black", face = "bold"),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(linewidth = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25)) +
  guides(color = guide_legend(keyheight = 1)) +
  labs(x = "Corrosion inhibitor type", color = "Inhibitor\nlevel",
       y = "Shannon alpha diversity score")
alpha.plot
ggsave("../Plots/alpha_inhibitorTypes.pdf", plot = alpha.plot, device = "pdf", width = 4, height = 4, units = "in")


