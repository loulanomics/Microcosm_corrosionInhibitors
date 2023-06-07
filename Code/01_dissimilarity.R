########################################
### Bray-Curtis dissimilarity analysis
### for Kimbell et al., 2023
### Lou LaMartina, finalized Jun 7, 2023
########################################


setwd("~/Desktop/16S/Kimbell/Final")
library(vegan)
library(pairwiseAdonis)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)


#################
### load data ###
#################
# removed A2_Day0_2 - low quality
# see 00_dada2.R for processing

# ASV counts; change hyphen to underscore
counts <- read.csv("Data/Kimbell2023_ASV_counts.csv", stringsAsFactors = F)
counts$Sample_name <- gsub("-", "_", counts$Sample_name)
rownames(counts) <- counts$Sample_name
counts <- counts[-1]


# sample metadata
info <- read.csv("Data/Kimbell2023_sample_info.csv", stringsAsFactors = F)
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

# stat
bray.matrix <- vegdist(relabun, "bray")
beta <- data.frame(ID_x = rownames(relabun), scores(bray.matrix))


# get comparisons
beta <- melt(beta, variable.name = "ID_y", value.name = "Bray", id.vars = "ID_x")
beta$ID_y <- as.character(beta$ID_y)
beta$Compare <- NA
for (i in 1:nrow(beta)) {
  beta$Compare[i] <- paste(sort(unlist(beta[i, c("ID_x", "ID_y")])),
                           collapse = "__")
}
beta <- unique(beta[beta$Bray > 0, c("Compare", "Bray")])


# extract info - x
beta$ID_x <- sapply(strsplit(beta$Compare, "__"), '[', 1)
beta$Set_x <- sapply(strsplit(beta$ID_x, "_"), '[', 1)

beta$Source_x <- "shallow"
beta$Source_x[beta$Set_x == "set1"] <- "deep"

beta$Inhib_x <- sapply(strsplit(beta$ID_x, "_"), '[', 2)
beta$Level_x <- "normal"
beta$Level_x[beta$Set_x == "set3"] <- "high"

beta$Day_x <- sapply(strsplit(beta$ID_x, "_"), '[', 3)
beta$Sample_x <- paste0(beta$Set_x, "_", beta$Inhib_x, "_", beta$Day_x)


# extract info - y
beta$ID_y <- sapply(strsplit(beta$Compare, "__"), '[', 2)
beta$Set_y <- sapply(strsplit(beta$ID_y, "_"), '[', 1)

beta$Source_y <- "shallow"
beta$Source_y[beta$Set_y == "set1"] <- "deep"

beta$Inhib_y <- sapply(strsplit(beta$ID_y, "_"), '[', 2)
beta$Level_y <- "normal"
beta$Level_y[beta$Set_y == "set3"] <- "high"

beta$Day_y <- sapply(strsplit(beta$ID_y, "_"), '[', 3)
beta$Sample_y <- paste0(beta$Set_y, "_", beta$Inhib_y, "_", beta$Day_y)
beta_all <- beta
#write.csv(beta_all, "Data/BrayCurtis_all.csv", row.names = F, na = "")




#####################
### means and SDs ###
#####################

# scientific notation after 10 decimal places
options(scipen = 10)


# compare within sets
beta <- beta_all[beta_all$Set_x == beta_all$Set_y,]


# first comparison is control
beta <- beta[beta$Inhib_x == "control",]


# compare same day
beta <- beta[beta$Day_x == beta$Day_y,]


# remove day zero
beta <- beta[beta$Day_x != "day0",]


# overall inhibs in shallow water (sets 2 and 3)
beta_shallow <- beta[beta$Set_x != "set1",]
beta_shallow <- data.frame(aggregate(Bray ~ Inhib_y, mean, data = beta_shallow),
                         aggregate(Bray ~ Inhib_y, sd, data = beta_shallow)[2])
colnames(beta_shallow) <- c("Inhibitor", "Mean", "SD")
beta_shallow$Set <- "set2+3"
beta_shallow$Level <- "normal+high"
beta_shallow$Source <- "shallow"


# inhibs at normal levels (sets 1 and 2)
beta_normal <- beta[beta$Set_x != "set3",]
beta_normal <- data.frame(aggregate(Bray ~ Inhib_y, mean, data = beta_normal),
                          aggregate(Bray ~ Inhib_y, sd, data = beta_normal)[2])
colnames(beta_normal) <- c("Inhibitor", "Mean", "SD")
beta_normal$Set <- "set1+2"
beta_normal$Level <- "normal"
beta_normal$Source <- "deep+shallow"


# beta for all
beta_summary <- data.frame(aggregate(Bray ~ Set_y + Source_y + Inhib_y + Level_y, mean, data = beta),
                           aggregate(Bray ~ Set_y + Source_y + Inhib_y + Level_y, sd, data = beta)[5])
colnames(beta_summary) <- c("Set", "Source", "Inhibitor", "Level", "Mean", "SD")


# combine
beta_summary <- bind_rows(beta_summary, beta_shallow, beta_normal)


# calculate CV
beta_summary$CV <- beta_summary$SD / beta_summary$Mean




#################
### normality ###
#################

beta_summary$Shapiro <- NA
for (s in unique(beta_summary$Set)) {
  for (i in unique(beta_summary$Inhibitor)) {
    
    if (s %in% c("set1", "set2", "set3")) {
      x <- beta[beta$Set_x == s & beta$Inhib_y == i,]
      x <- shapiro.test(x$Bray)$p.value
      beta_summary$Shapiro[beta_summary$Set == s & beta_summary$Inhibitor == i] <- x
    }
    
    if (s == "set2+3") {
      x <- beta[beta$Set_x != "set1" & beta$Inhib_y == i,]
      x <- shapiro.test(x$Bray)$p.value
      beta_summary$Shapiro[beta_summary$Set == s & beta_summary$Inhibitor == i] <- x
    }
    
    if (s == "set1+2") {
      x <- beta[beta$Set_x != "set3" & beta$Inhib_y == i,]
      x <- shapiro.test(x$Bray)$p.value
      beta_summary$Shapiro[beta_summary$Set == s & beta_summary$Inhibitor == i] <- x
    }
  }
}



# two sig figs
beta_summary$Mean <- signif(beta_summary$Mean, 2)
beta_summary$SD <- signif(beta_summary$SD, 2)
beta_summary$CV <- signif(beta_summary$CV, 2)
beta_summary$Shapiro <- signif(beta_summary$Shapiro, 2)
#write.csv(beta_summary, "Data/BrayCurtis_summary.csv", row.names = F, na = "")



###################
### comparisons ###
###################

options(scipen = 0)
inhibs <- c("Na2SiO3", "NaPO4", "ZnPO4")


# inhibs overall
for (i in inhibs) {
  cat("\n", i, "- ")
  cat(signif(wilcox.test(subset(beta, Inhib_y == i)$Bray,
                         subset(beta, Inhib_y == "control")$Bray,
                         paired = F, alternative = "greater")$p.value, 2))
}
# Na2SiO3 - 4.7e-05
# NaPO4 - 1.9e-06
# ZnPO4 - 2.7e-10


# inhib levels overall
for (i in c("normal", "high")) {
  cat("\n", i, "- ")
  cat(signif(wilcox.test(subset(beta, Level_x == i & Inhib_y != "control")$Bray,
                         subset(beta, Level_x == i & Inhib_y == "control")$Bray,
                         paired = F, alternative = "greater")$p.value, 2))
}
# normal - 2.3e-07
# high - 0.0028


# inhibs in set 3
for (i in inhibs) {
  cat("\n", i, "- ")
  cat(signif(wilcox.test(subset(beta, Set_x == "set3" & Inhib_y == i)$Bray,
                           subset(beta, Set_x == "set3" & Inhib_y == "control")$Bray,
                           paired = F, alternative = "greater")$p.value, 2))
}
# Na2SiO3 - 7.4e-06
# NaPO4 - 0.39
# ZnPO4 - 7.4e-06




#################
### PERMANOVA ###
#################

# permanova - controlling for water source
do.call(rbind, pairwise.adonis2(relabun ~ Inhibitor, strata = "Water_source", data = info)[-1])
#                              Df SumOfSqs      R2      F Pr(>F)    
# control_vs_Na2SiO3.Inhibitor  1   0.4744 0.04159 1.8226  0.005 ** 
# control_vs_Na2SiO3.Residual  42  10.9334 0.95841                  
# control_vs_Na2SiO3.Total     43  11.4078 1.00000                  
# control_vs_NaPO4.Inhibitor    1   0.3023 0.02747 1.1863  0.053 .  
# control_vs_NaPO4.Residual    42  10.7014 0.97253                  
# control_vs_NaPO4.Total       43  11.0037 1.00000                  
# control_vs_ZnPO4.Inhibitor    1   0.9284 0.07577 3.4432  0.001 ***
# control_vs_ZnPO4.Residual    42  11.3247 0.92423                  
# control_vs_ZnPO4.Total       43  12.2532 1.00000                  
# Na2SiO3_vs_NaPO4.Inhibitor    1   0.4277 0.04474 1.5923  0.019 *  
# Na2SiO3_vs_NaPO4.Residual    34   9.1319 0.95526                  
# Na2SiO3_vs_NaPO4.Total       35   9.5595 1.00000                  
# Na2SiO3_vs_ZnPO4.Inhibitor    1   0.6560 0.06301 2.2862  0.002 ** 
# Na2SiO3_vs_ZnPO4.Residual    34   9.7552 0.93699                  
# Na2SiO3_vs_ZnPO4.Total       35  10.4111 1.00000                  
# NaPO4_vs_ZnPO4.Inhibitor      1   0.7132 0.06967 2.5462  0.001 ***
# NaPO4_vs_ZnPO4.Residual      34   9.5232 0.93033                  
# NaPO4_vs_ZnPO4.Total         35  10.2364 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# permanova - controlling for water source
perman <- adonis2(bray.matrix ~ Inhibitor * Strength + Replicate + Day, 
                  strata = info$Water_source, data = info)


# extract info
perman <- data.frame(Variables = rownames(perman), perman)
perman[-1] <- apply(perman[-1], 2, function(x) signif(x, 2))


# post-hoc pairwise multilevel comparison
posthoc <- pairwise.adonis2(relabun ~ Inhibitor * Strength * Replicate * Day, 
                            strata = "Water_source", data = info)
posthoc <- do.call(rbind, posthoc[-1])
posthoc <- data.frame(Comparison = sapply(strsplit(rownames(posthoc), "\\."), '[', 1),
                      Variables = sapply(strsplit(rownames(posthoc), "\\."), '[', 2),
                      posthoc)


# keep p values
posthoc$Pr..F. <- signif(posthoc$Pr..F., 2)
posthoc.p <- dcast(Variables ~ Comparison, value.var = "Pr..F.", data = posthoc)
colnames(posthoc.p)[-1] <- paste0(colnames(posthoc.p)[-1], ".p")


# keep R squared values
posthoc$R2 <- signif(posthoc$R2, 2)
posthoc.R2 <- dcast(Variables ~ Comparison, value.var = "R2", data = posthoc)
colnames(posthoc.R2)[-1] <- paste0(colnames(posthoc.p)[-1], ".R2")


# combine p and R2
posthoc.p <- merge(posthoc.p, posthoc.R2, by = "Variables")
posthoc.p <- posthoc.p[c("Variables", sort(colnames(posthoc.p)[-1]))]


# combine permanova
perman <- merge(perman, posthoc.p, by = "Variables")


# validate with analysis of multivariate homogeneity of group dispersions (variances)
inhibdisp <- betadisper(bray.matrix, info$Inhibitor) 
anova(inhibdisp) # p = 0.01875 *
leveldisp <- betadisper(bray.matrix, info$Strength) 
anova(leveldisp) # p = 0.03213 *


# confidence intervals on the differences between the mean distance-to-centroid
plot(TukeyHSD(inhibdisp))




##################
### ordination ###
##################

# stat & extract info
nmds <- metaMDS(relabun, distance = "bray")
nmds$stress # 0.08367593
nmds <- data.frame(Sample_ID = rownames(scores(nmds)$sites), scores(nmds)$sites)
nmds <- merge(nmds, info, by = "Sample_ID")
#write.csv(nmds, "Data/metaMDS.csv", row.names = F, na = "")


# labels
inhiblabs <- c(Na2SiO3 = expression(Na[2]*SiO[3]),
               NaPO4 = expression(NaPO[4]),
               ZnPO4 = expression(ZnPO[4]))


# plot
nmds.plot <-
  ggplot(nmds, aes(x = NMDS1, y = NMDS2, color = Inhibitor, shape = Set)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.2) +
  geom_point(size = 2) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  scale_shape_manual(values = c(17, 19, 15),
                     labels = c(set1 = "deep + normal", 
                                set2 = "shallow + normal", 
                                set3 = "shallow + high")) +
  scale_color_manual(values = c("grey70", brewer.pal(3, "Dark2")),
                     labels = inhiblabs) +
  annotate("text", x = -1.6, y = 1.5, label = "stress = 0.084", size = 2.75, color = "grey60") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        legend.title = element_text(size = 10, color = "black", face = "bold"),
        legend.text = element_text(size = 9, color = "black", hjust = 0),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.25)) +
  guides(color = guide_legend(keyheight = 0.75), shape = guide_legend(keyheight = 0.75)) +
  labs(shape = "Water source +\ninhibitor level")
nmds.plot

#ggsave("PLOTS/nmds.pdf", plot = nmds.plot, device = "pdf", width = 6, height = 3, units = "in")




############
### plot ###
############

# effects on inhibitor levels
beta_levels <- beta[beta$Set_x != "set1", 
                    c("Bray", "Set_y", "Source_y", "Inhib_y", "Level_y", "Day_y")]
colnames(beta_levels)[-1] <- c("Set", "Source", "Inhibitor", "Level", "Day")
levels_summary <- beta_summary[beta_summary$Set %in% c("set2", "set3"),]


# plot
levels.plot <-
  ggplot(beta_levels, aes(x = Inhibitor, y = Bray, color = Level)) +
  
  # points
  geom_jitter(shape = 1, size = 0.5, position = 
                position_jitterdodge(jitter.width = 0.1, dodge.width = 0.25)) +
  geom_errorbar(data = levels_summary, show.legend = F,
                aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD), 
                width = 0.2, linewidth = 1, position = position_dodge(0.25)) +
  geom_point(data = levels_summary, aes(y = Mean, shape = Level), size = 3, 
             position = position_dodge(0.25)) +
  
  scale_x_discrete(labels = inhiblabs) +
  scale_color_manual(values = c("#c02f1d", "#0d3356")) +
  scale_shape_manual(values = c(17, 19)) +
  scale_y_continuous(limits = c(0.14,0.78), expand = c(0,0.02)) +
  
  # theme
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 9, color = "black", face = "bold"),
        axis.title.y = element_text(size = 9, color = "black", face = "bold"),
        legend.position = "top",
        legend.margin = margin(b = -0.3, unit = "cm"),
        legend.title = element_text(size = 9, color = "black", face = "bold"),
        strip.text = element_text(size = 8, color = "black", face = "bold"),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(linewidth = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25)) +
  guides(color = guide_legend(title.hjust = 0.5, keyheight = 0.75, nrow = 2, title.position = "left")) +
  labs(x = "Corrosion inhibitor", color = "Inhibitor level", shape = "Inhibitor level",
       y = "Bray-Curtis dissimilarity score\ncompared to same-day control")
levels.plot

#ggsave("PLOTS/beta_inhibitorLevels.pdf", plot = levels.plot, device = "pdf", width = 3, height = 4, units = "in")


# effects on deep water
beta_deep <- beta[beta$Set_x != "set3", 
                   c("Bray", "Set_y", "Source_y", "Inhib_y", "Level_y", "Day_y")]
colnames(beta_deep)[-1] <- c("Set", "Source", "Inhibitor", "Level", "Day")
treat_summary <- beta_summary[beta_summary$Set %in% c("set1", "set2"),]


# plot
deep.plot <-
  ggplot(beta_deep, aes(x = Inhibitor, y = Bray, color = Source)) +
  
  # points
  geom_jitter(shape = 1, size = 0.5, position = 
                position_jitterdodge(jitter.width = 0.1, dodge.width = 0.25)) +
  geom_errorbar(data = treat_summary, show.legend = F,
                aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD), 
                width = 0.2, linewidth = 1, position = position_dodge(0.25)) +
  geom_point(data = treat_summary, aes(y = Mean, shape = Source), size = 3, 
             position = position_dodge(0.25)) +
  
  scale_x_discrete(labels = inhiblabs) +
  scale_color_manual(values = c("#d3b53d", "#0d3356")) +
  scale_shape_manual(values = c(15, 19)) +
  scale_y_continuous(limits = c(0.14,0.78), expand = c(0,0.02)) +
  
  # theme
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 9, color = "black", face = "bold"),
        axis.title.y = element_text(size = 9, color = "black", face = "bold"),
        legend.position = "top",
        legend.margin = margin(b = -0.3, unit = "cm"),
        legend.title = element_text(size = 9, color = "black", face = "bold"),
        strip.text = element_text(size = 8, color = "black", face = "bold"),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(linewidth = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25)) +
  guides(color = guide_legend(title.hjust = 0.5, keyheight = 0.75, nrow = 2, title.position = "right")) +
  labs(x = "Corrosion inhibitor", color = "Water source", shape = "Water source",
       y = "Bray-Curtis dissimilarity score\ncompared to same-day control")
deep.plot

#ggsave("PLOTS/beta_inhibitorTreatment.pdf", plot = deep.plot, device = "pdf", width = 3, height = 4, units = "in")




############
### save ###
############

# compare within same sets
#write.csv(beta, "Data/BrayCurtis.csv", row.names = F, na = "")


# permanova & posthoc results
#write.csv(perman, "Data/permanova.csv", row.names = F, na = "")
