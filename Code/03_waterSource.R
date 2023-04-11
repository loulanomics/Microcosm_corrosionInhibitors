#########################################
### effects of corrosion inhibitor types
### on microbial communities from treated
### and untreated water sources
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

# normal inhibitor levels, treated (set1) vs untreated (set2) source water
# compared same day, same set, to same day & set controls
beta <- rbind(
  beta[beta$Set1 == "set1" & beta$Set2 == "set1" &
         beta$Treat1 == "control" & 
         beta$Day1 == beta$Day2, -1],
  beta[beta$Set1 == "set2" & beta$Set2 == "set2" &
         beta$Treat1 == "control" & 
         beta$Day1 == beta$Day2, -1])


# don't want day zero; comparing treatments to controls on same day
beta <- beta[beta$Day1 != "day0",]
#beta$Treat3 <- beta$Treat2
#beta$Treat3[beta$Treat3 != "control"] <- "treatment"


# summary stats for plot
beta.stats <- data.frame(aggregate(Bray ~ Set1 + Treat2, mean, data = beta),
                         aggregate(Bray ~ Set1 + Treat2, sd, data = beta)[3])
colnames(beta.stats)[3:4] <- c("mean", "sd")




###############
### overall ###
###############

###################
### normality tests

### treated

# beta-curtis btwn replicates of same-day controls
treatctrl.beta <- beta[beta$Treat2 == "control" & beta$Set1 == "set1",]
cat(round(mean(treatctrl.beta$Bray), digits = 2), "±", round(sd(treatctrl.beta$Bray), digits = 2))
shapiro.verbose(treatctrl.beta$Bray)
# 0.16 ± 0.01
# p = 0.2127929 - not normal distribution - use wilcox.test


# beta-curtis btwn replicates of same-day controls
treatexp.beta <- beta[beta$Treat2 != "control" & beta$Set1 == "set1",]
cat(round(mean(treatexp.beta$Bray), digits = 2), "±", round(sd(treatexp.beta$Bray), digits = 2))
shapiro.verbose(treatexp.beta$Bray)
# 0.29 ± 0.06
# p = 0.1774276 - not normal distribution - use wilcox.test


# beta-curtis btwn Na2SiO3 samples and same-day controls
treatss.beta <- beta[beta$Treat2 == "Na2SiO3" & beta$Set1 == "set1",]
cat(round(mean(treatss.beta$Bray), digits = 2), "±", round(sd(treatss.beta$Bray), digits = 2))
shapiro.verbose(treatss.beta$Bray)
# 0.24 ± 0.04
# p = 0.7301303 - not normal distribution - use wilcox.test


# beta-curtis btwn Na2SiO3 samples and same-day controls
treatso.beta <- beta[beta$Treat2 == "NaPO4" & beta$Set1 == "set1",]
cat(round(mean(treatso.beta$Bray), digits = 2), "±", round(sd(treatso.beta$Bray), digits = 2))
shapiro.verbose(treatso.beta$Bray)
# 0.28 ± 0.04
# p = 0.6429402 - not normal distribution - use wilcox.test


# beta-curtis btwn Na2SiO3 samples and same-day controls
treatzo.beta <- beta[beta$Treat2 == "ZnPO4" & beta$Set1 == "set1",]
cat(round(mean(treatzo.beta$Bray), digits = 2), "±", round(sd(treatzo.beta$Bray), digits = 2))
shapiro.verbose(treatzo.beta$Bray)
# 0.36 ± 0.04
# p = 0.9150062 - not normal distribution - use wilcox.test


### untreated

# beta-curtis btwn replicates of same-day controls
untreatctrl.beta <- beta[beta$Treat2 == "control" & beta$Set1 == "set2",]
cat(round(mean(untreatctrl.beta$Bray), digits = 2), "±", round(sd(untreatctrl.beta$Bray), digits = 2))
shapiro.verbose(untreatctrl.beta$Bray)
# 0.21 ± 0.04
# p = 0.5197797 - not normal distribution - use wilcox.test


# beta-curtis btwn replicates of same-day controls
untreatexp.beta <- beta[beta$Treat2 != "control" & beta$Set1 == "set2",]
cat(round(mean(untreatexp.beta$Bray), digits = 2), "±", round(sd(untreatexp.beta$Bray), digits = 2))
shapiro.verbose(untreatexp.beta$Bray)
# 0.29 ± 0.06
# p = 0.1774276 - not normal distribution - use wilcox.test


# beta-curtis btwn Na2SiO3 samples and same-day controls
untreatss.beta <- beta[beta$Treat2 == "Na2SiO3" & beta$Set1 == "set2",]
cat(round(mean(untreatss.beta$Bray), digits = 2), "±", round(sd(untreatss.beta$Bray), digits = 2))
shapiro.verbose(untreatss.beta$Bray)
# 0.26 ± 0.04
# p = 0.8423096 - not normal distribution - use wilcox.test


# beta-curtis btwn NaPO4 samples and same-day controls
untreatso.beta <- beta[beta$Treat2 == "NaPO4" & beta$Set1 == "set2",]
cat(round(mean(untreatso.beta$Bray), digits = 2), "±", round(sd(untreatso.beta$Bray), digits = 2))
shapiro.verbose(untreatso.beta$Bray)
# 0.36 ± 0.06
# p = 0.822262 - not normal distribution - use wilcox.test


# beta-curtis btwn ZnPO4 samples and same-day controls
untreatzo.beta <- beta[beta$Treat2 == "ZnPO4" & beta$Set1 == "set2",]
cat(round(mean(untreatzo.beta$Bray), digits = 2), "±", round(sd(untreatzo.beta$Bray), digits = 2))
shapiro.verbose(untreatzo.beta$Bray)
# 0.63 ± 0.02
# p = 0.7176429 - not normal distribution - use wilcox.test



###############
### comparisons

### treated water

# were exp samples a greater distance from same-day controls
# than controls to replicate controls?
cat("p =", wilcox.test(x = treatexp.beta$Bray, y = treatctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# 3.829394e-05


# were distances btwn treated water Na2SiO3-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = treatss.beta$Bray, y = treatctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 0.01121876


# were distances btwn treated water NaPO4-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = treatso.beta$Bray, y = treatctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 2.971857e-05


# were distances btwn treated water ZnPO4-amended samples and same-day controls
# greater than replicate controls to each other?
cat("p =", wilcox.test(x = treatzo.beta$Bray, y = treatctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# p = 7.429641e-06


### untreated water

# were exp samples a greater distance from same-day controls
# than controls to replicate controls?
cat("p =", wilcox.test(x = untreatexp.beta$Bray, y = untreatctrl.beta$Bray, paired = F, alternative = "greater")$p.value)
# 0.0002936866


# ZnPO4 more of an impact in untreaeted than treated?
cat("p =", wilcox.test(x = untreatzo.beta$Bray, y = treatzo.beta$Bray, paired = F, alternative = "greater")$p.value)
# 1.101912e-10


# NaPO4 more of an impact in untreaeted than treated?
cat("p =", wilcox.test(x = untreatso.beta$Bray, y = treatso.beta$Bray, paired = F, alternative = "greater")$p.value)
# 7.738232e-05


# Na2SiO3 more of an impact in untreaeted than treated?
cat("p =", wilcox.test(x = untreatss.beta$Bray, y = treatss.beta$Bray, paired = F, alternative = "greater")$p.value)
# 0.04535378



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
  scale_color_manual(labels = c(set1 = "treated", set2 = "untreated"),
                     values = brewer.pal(8, "Set1")[3:2]) +
  
  # theme
  theme_classic() +
  theme(axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        legend.position = c(0.17,0.85),
        legend.background = element_rect(fill = NA, color = "grey80"),
        legend.title = element_text(size = 9, color = "black", face = "bold"),
        strip.text = element_text(size = 10, color = "black", face = "bold"),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(linewidth = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25)) +
  guides(color = guide_legend(keyheight = 1)) +
  labs(x = "Corrosion inhibitor type", color = "Water\nsource",
       y = "Bray-Curtis dissimilarity score\ncompared to same-day control")
beta.plot
ggsave("../Plots/beta_waterSource.pdf", plot = beta.plot, device = "pdf", width = 4, height = 4, units = "in")




#######################
### alpha diversity ###
#######################

########
### stat 

# stat
alpha <- data.frame(Sample_ID = rownames(relabun), Shannon = diversity(relabun, "shannon"))


# add info
alpha <- merge(alpha, info, by = "Sample_ID")
alpha <- alpha[alpha$Set != "set3" & alpha$Day != "day0",]


# summary stats for plot
alpha.stats <- data.frame(aggregate(Shannon ~ Set + Inhibitor, mean, data = alpha),
                          aggregate(Shannon ~ Set + Inhibitor, sd, data = alpha)[3])
colnames(alpha.stats)[3:4] <- c("mean", "sd")



###################
### normality tests

### treated water

# shannon diversity in controls
treatctrl.alpha <- alpha[alpha$Set == "set1" & alpha$Treatment == "control",]
cat(round(mean(treatctrl.alpha$Shannon), digits = 2), "±", round(sd(treatctrl.alpha$Shannon), digits = 2))
shapiro.verbose(treatctrl.alpha$Shannon) 
# 5.52 ± 0.24
# p = 0.1868189 - not normal distribution - use wilcox.test


# shannon diversity in Na2SiO3
treatss.alpha <- alpha[alpha$Set == "set1" & alpha$Inhibitor == "Na2SiO3",]
cat(round(mean(treatss.alpha$Shannon), digits = 2), "±", round(sd(treatss.alpha$Shannon), digits = 2))
shapiro.verbose(treatss.alpha$Shannon) 
# 5.62 ± 0.14
# p = 0.2082664 - not normal distribution - use wilcox.test


# shannon diversity in NaPO4
treatso.alpha <- alpha[alpha$Set == "set1" & alpha$Inhibitor == "NaPO4",]
cat(round(mean(treatso.alpha$Shannon), digits = 2), "±", round(sd(treatso.alpha$Shannon), digits = 2))
shapiro.verbose(treatso.alpha$Shannon) 
# 5.67 ± 0.06
# p = 0.1594336 - not normal distribution - use wilcox.test


# shannon diversity in ZnPO4
treatzo.alpha <- alpha[alpha$Set == "set1" & alpha$Inhibitor == "ZnPO4",]
cat(round(mean(treatzo.alpha$Shannon), digits = 2), "±", round(sd(treatzo.alpha$Shannon), digits = 2))
shapiro.verbose(treatzo.alpha$Shannon) 
# 5.45 ± 0.16
# p = 0.2905456 - not normal distribution - use wilcox.test



### untreated water

# shannon diversity in controls
untreatctrl.alpha <- alpha[alpha$Set == "set2" & alpha$Treatment == "control",]
cat(round(mean(untreatctrl.alpha$Shannon), digits = 2), "±", round(sd(untreatctrl.alpha$Shannon), digits = 2))
shapiro.verbose(untreatctrl.alpha$Shannon) 
# 5.54 ± 0.21
# p = 0.305246 - not normal distribution - use wilcox.test


# shannon diversity in Na2SiO3
untreatss.alpha <- alpha[alpha$Set == "set2" & alpha$Inhibitor == "Na2SiO3",]
cat(round(mean(untreatss.alpha$Shannon), digits = 2), "±", round(sd(untreatss.alpha$Shannon), digits = 2))
shapiro.verbose(untreatss.alpha$Shannon) 
# 5.52 ± 0.34
# p = 0.3759599 - not normal distribution - use wilcox.test


# shannon diversity in NaPO4
untreatso.alpha <- alpha[alpha$Set == "set2" & alpha$Inhibitor == "NaPO4",]
cat(round(mean(untreatso.alpha$Shannon), digits = 2), "±", round(sd(untreatso.alpha$Shannon), digits = 2))
shapiro.verbose(untreatso.alpha$Shannon) 
# 5.66 ± 0.18
# p = 0.07938339 - not normal distribution - use wilcox.test


# shannon diversity in ZnPO4
untreatzo.alpha <- alpha[alpha$Set == "set2" & alpha$Inhibitor == "ZnPO4",]
cat(round(mean(untreatzo.alpha$Shannon), digits = 2), "±", round(sd(untreatzo.alpha$Shannon), digits = 2))
shapiro.verbose(untreatzo.alpha$Shannon) 
# 5.2 ± 0.23
# p = 0.08383906 - not normal distribution - use wilcox.test



###############
### comparisons

### treated water

# differences in community richness between Na2SiO3 samples and controls?
cat("p =", wilcox.test(x = treatctrl.alpha$Shannon, y = treatss.alpha$Shannon, paired = F)$p.value)
# p = 0.4848485


# differences in community richness between NaPO4 samples and controls?
cat("p =", wilcox.test(x = treatctrl.alpha$Shannon, y = treatso.alpha$Shannon, paired = F)$p.value)
# p = 0.4848485


# differences in community richness between ZnPO4 samples and controls?
cat("p =", wilcox.test(x = treatctrl.alpha$Shannon, y = treatzo.alpha$Shannon, paired = F)$p.value)
# p = 0.6991342



### untreated

# differences in community richness between Na2SiO3 samples and controls?
cat("p =", wilcox.test(x = untreatctrl.alpha$Shannon, y = untreatss.alpha$Shannon, paired = F)$p.value)
# p = 0.9372294


# differences in community richness between NaPO4 samples and controls?
cat("p =", wilcox.test(x = untreatctrl.alpha$Shannon, y = untreatso.alpha$Shannon, paired = F)$p.value)
# p = 0.2402597


# normal inhibitor levels - 
# differences in community richness between ZnPO4 samples and controls?
cat("p =", wilcox.test(x = untreatctrl.alpha$Shannon, y = untreatzo.alpha$Shannon, paired = F, alternative = "greater")$p.value)
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
  scale_color_manual(labels = c(set1 = "treated", set2 = "untreated"),
                     values = brewer.pal(8, "Set1")[3:2]) +
  
  # theme
  theme_classic() +
  theme(axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        legend.position = c(0.17,0.15),
        legend.background = element_rect(fill = NA, color = "grey80"),
        legend.title = element_text(size = 9, color = "black", face = "bold"),
        strip.text = element_text(size = 10, color = "black", face = "bold"),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(linewidth = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25)) +
  guides(color = guide_legend(keyheight = 1)) +
  labs(x = "Corrosion inhibitor type", color = "Water\nsource",
       y = "Shannon alpha diversity score")
alpha.plot
ggsave("../Plots/alpha_waterSource.pdf", plot = alpha.plot, device = "pdf", width = 4, height = 4, units = "in")
