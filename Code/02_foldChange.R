########################################
### ASV fold-change to treatments
### for Kimbell et al., 2023
### Lou LaMartina, finalized Jun 7, 2023
########################################


setwd("~/Desktop/16S/Kimbell/Final")
library(ggplot2)
library(reshape2)



#################
### load data ###
#################

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


# taxa
taxa <- read.csv("Data/Kimbell2023_ASV_taxa.csv", stringsAsFactors = F)
taxa[taxa == ""] <- NA


# only ASVs in all replicates
reps <- list()
for (i in unique(info$Sample_name)) {
  df <- counts[rownames(counts) %in% info$Sample_ID[info$Sample_name == i],]
  df[df > 0] <- 1
  reps[[i]] <- names(which(colSums(df) == nrow(df)))
}; rm(df)
reps <- sort(unique(unlist(reps)))


# convert to relative abundance
reps <- counts[reps]
relabun <-  reps / rowSums(reps)
rm(reps)
ncol(relabun) # 2728


# get min relative abundance
minval <- min(melt(relabun)$value[melt(relabun)$value != 0])




##################
### label taxa ###
##################

# make labels
labs <- taxa


# remove non specifics
for (i in colnames(labs)) {
  labs[[i]][grep("group", labs[[i]], ignore.case = T)] <- NA
  labs[[i]][grep("^clade", labs[[i]], ignore.case = T)] <- NA
  labs[[i]][grep("endosymbionts", labs[[i]], ignore.case = T)] <- NA
  labs[[i]][grep("unknown", labs[[i]], ignore.case = T)] <- NA
  labs[[i]][grep("possible", labs[[i]], ignore.case = T)] <- NA
  labs[[i]][grep("family", labs[[i]], ignore.case = T)] <- NA
  labs[[i]][grep("rhizobium-", labs[[i]], ignore.case = T)] <- "Rhizobium"
  labs[[i]][grep("-", labs[[i]])] <- NA
  labs[[i]][grep("\\.", labs[[i]])] <- NA
  labs[[i]] <- gsub("Candidatus_", "", labs[[i]], ignore.case = F)
  labs[[i]] <- gsub("_clade", " clade", labs[[i]], ignore.case = F)
  labs[[i]] <- gsub("_lineage", " lineage", labs[[i]], ignore.case = F)
  labs[[i]] <- gsub("_sensu.*", "", labs[[i]], ignore.case = F)
  labs[[i]] <- gsub("_Incertae.*", "", labs[[i]], ignore.case = F)
  labs[[i]] <- gsub("_[[:digit:]].*", "", labs[[i]], ignore.case = F)
}

rmv <- unique(unlist(labs[3:7])[grep("[[:digit:]]", unlist(labs[3:7]))])
rmv <- rmv[grep(" ", rmv, invert = T)]

for (i in colnames(labs)) {
  labs[[i]][labs[[i]] %in% rmv] <- NA
}


# individual changes
labs$Phylum[labs$Class == "Deinococci"] <- "Deinococcota"
labs$Phylum[is.na(labs$Phylum)] <- labs$Kingdom[is.na(labs$Phylum)]
labs$Class[is.na(labs$Class) == F] <- paste0("c. ", labs$Class[is.na(labs$Class) == F])
labs$Order[is.na(labs$Order) == F] <- paste0("o. ", labs$Order[is.na(labs$Order) == F])
labs$Family[is.na(labs$Family) == F] <- paste0("f. ", labs$Family[is.na(labs$Family) == F])
labs$Genus[is.na(labs$Genus) == F] <- paste0("g. ", labs$Genus[is.na(labs$Genus) == F])


# format
labs$Label <- NA
for (i in 1:nrow(labs)) {
  lab <- unlist(labs[i, 4:7])
  lab <- lab[is.na(lab) == F]
  lab <- lab[length(lab)]
  labs$Label[i] <- paste0(labs$Phylum[i], " (", lab, ")")
}

labs$Label <- gsub(" \\()", "", labs$Label)
taxa$Label <- labs$Label


# update file
#write.csv(taxa, "Data/Kimbell2023_ASV_taxa.csv", row.names = F, na = "")




#############
### set 1 ###
#############

# subset set1, day 7
set1_info <- subset(info, Set == "set1" & Day == "day7")
set1_relabun <- subset(relabun, rownames(relabun) %in% set1_info$Sample_ID)
set1_relabun <- set1_relabun[colSums(set1_relabun) > 0]
ncol(set1_relabun) # 1378
set1_relabun[set1_relabun == 0] <- minval


# controls
set1_ctrl_relabun <- set1_relabun[grep("control", rownames(set1_relabun)),]
set1_ctrl.m <- melt(data.frame(SampleID = rownames(set1_ctrl_relabun), set1_ctrl_relabun),
                    value.name = "Ctrl", variable.name = "ASV")
set1_ctrl.m$Rep <- sapply(strsplit(set1_ctrl.m$SampleID, "_"), '[', 4)
set1_ctrl.m$Set <- "set1"
set1_ctrl.m <- set1_ctrl.m[c("Set", "Rep", "ASV", "Ctrl")]


# Na2SiO3
set1_ss_relabun <- set1_relabun[grep("Na2SiO3", rownames(set1_relabun)),]
set1_ss.m <- melt(data.frame(SampleID = rownames(set1_ss_relabun), set1_ss_relabun),
                  value.name = "Exp", variable.name = "ASV")
set1_ss.m$Rep <- sapply(strsplit(set1_ss.m$SampleID, "_"), '[', 4)
set1_ss.m$Inhibitor <- "Na2SiO3"
set1_ss.m$Set <- "set1"
set1_ss.m <- set1_ss.m[c("Set", "Rep", "Inhibitor", "ASV", "Exp")]
identical(set1_ss.m$ASV, set1_ctrl.m$ASV); identical(set1_ss.m$Rep, set1_ctrl.m$Rep)
set1_ss.m$Ctrl <- set1_ctrl.m$Ctrl


# NaPO4
set1_so_relabun <- set1_relabun[grep("NaPO4", rownames(set1_relabun)),]
set1_so.m <- melt(data.frame(SampleID = rownames(set1_so_relabun), set1_so_relabun),
                  value.name = "Exp", variable.name = "ASV")
set1_so.m$Rep <- sapply(strsplit(set1_so.m$SampleID, "_"), '[', 4)
set1_so.m$Inhibitor <- "NaPO4"
set1_so.m$Set <- "set1"
set1_so.m <- set1_so.m[c("Set", "Rep", "Inhibitor", "ASV", "Exp")]
identical(set1_so.m$ASV, set1_ctrl.m$ASV); identical(set1_so.m$Rep, set1_ctrl.m$Rep)
set1_so.m$Ctrl <- set1_ctrl.m$Ctrl


# ZnPO4
set1_zo_relabun <- set1_relabun[grep("ZnPO4", rownames(set1_relabun)),]
set1_zo.m <- melt(data.frame(SampleID = rownames(set1_zo_relabun), set1_zo_relabun),
                  value.name = "Exp", variable.name = "ASV")
set1_zo.m$Rep <- sapply(strsplit(set1_zo.m$SampleID, "_"), '[', 4)
set1_zo.m$Inhibitor <- "ZnPO4"
set1_zo.m$Set <- "set1"
set1_zo.m <- set1_zo.m[c("Set", "Rep", "Inhibitor", "ASV", "Exp")]
identical(set1_zo.m$ASV, set1_ctrl.m$ASV); identical(set1_zo.m$Rep, set1_ctrl.m$Rep)
set1_zo.m$Ctrl <- set1_ctrl.m$Ctrl


# combine
set1_fold <- rbind(set1_ss.m, set1_so.m, set1_zo.m)


# fold change
set1_fold$Diff[set1_fold$Ctrl < set1_fold$Exp] <- "increase"
set1_fold$Diff[set1_fold$Ctrl > set1_fold$Exp] <- "decrease"
set1_fold$logCtrl <- log2(set1_fold$Ctrl)
set1_fold$logExp <- log2(set1_fold$Exp)
set1_fold$Fold <- set1_fold$logExp - set1_fold$logCtrl




#############
### set 2 ###
#############

# subset set2, day 7
set2_info <- subset(info, Set == "set2" & Day == "day7")
set2_relabun <- subset(relabun, rownames(relabun) %in% set2_info$Sample_ID)
set2_relabun <- set2_relabun[colSums(set2_relabun) > 0]
ncol(set2_relabun) # 1447
set2_relabun[set2_relabun == 0] <- minval


# controls
set2_ctrl_relabun <- set2_relabun[grep("control", rownames(set2_relabun)),]
set2_ctrl.m <- melt(data.frame(SampleID = rownames(set2_ctrl_relabun), set2_ctrl_relabun),
                    value.name = "Ctrl", variable.name = "ASV")
set2_ctrl.m$Rep <- sapply(strsplit(set2_ctrl.m$SampleID, "_"), '[', 4)
set2_ctrl.m$Set <- "set2"
set2_ctrl.m <- set2_ctrl.m[c("Set", "Rep", "ASV", "Ctrl")]


# Na2SiO3
set2_ss_relabun <- set2_relabun[grep("Na2SiO3", rownames(set2_relabun)),]
set2_ss.m <- melt(data.frame(SampleID = rownames(set2_ss_relabun), set2_ss_relabun),
                  value.name = "Exp", variable.name = "ASV")
set2_ss.m$Rep <- sapply(strsplit(set2_ss.m$SampleID, "_"), '[', 4)
set2_ss.m$Inhibitor <- "Na2SiO3"
set2_ss.m$Set <- "set2"
set2_ss.m <- set2_ss.m[c("Set", "Rep", "Inhibitor", "ASV", "Exp")]
identical(set2_ss.m$ASV, set2_ctrl.m$ASV); identical(set2_ss.m$Rep, set2_ctrl.m$Rep)
set2_ss.m$Ctrl <- set2_ctrl.m$Ctrl


# NaPO4
set2_so_relabun <- set2_relabun[grep("NaPO4", rownames(set2_relabun)),]
set2_so.m <- melt(data.frame(SampleID = rownames(set2_so_relabun), set2_so_relabun),
                  value.name = "Exp", variable.name = "ASV")
set2_so.m$Rep <- sapply(strsplit(set2_so.m$SampleID, "_"), '[', 4)
set2_so.m$Inhibitor <- "NaPO4"
set2_so.m$Set <- "set2"
set2_so.m <- set2_so.m[c("Set", "Rep", "Inhibitor", "ASV", "Exp")]
identical(set2_so.m$ASV, set2_ctrl.m$ASV); identical(set2_so.m$Rep, set2_ctrl.m$Rep)
set2_so.m$Ctrl <- set2_ctrl.m$Ctrl


# ZnPO4
set2_zo_relabun <- set2_relabun[grep("ZnPO4", rownames(set2_relabun)),]
set2_zo.m <- melt(data.frame(SampleID = rownames(set2_zo_relabun), set2_zo_relabun),
                  value.name = "Exp", variable.name = "ASV")
set2_zo.m$Rep <- sapply(strsplit(set2_zo.m$SampleID, "_"), '[', 4)
set2_zo.m$Inhibitor <- "ZnPO4"
set2_zo.m$Set <- "set2"
set2_zo.m <- set2_zo.m[c("Set", "Rep", "Inhibitor", "ASV", "Exp")]
identical(set2_zo.m$ASV, set2_ctrl.m$ASV); identical(set2_zo.m$Rep, set2_ctrl.m$Rep)
set2_zo.m$Ctrl <- set2_ctrl.m$Ctrl


# combine
set2_fold <- rbind(set2_ss.m, set2_so.m, set2_zo.m)


# fold change
set2_fold$Diff[set2_fold$Ctrl < set2_fold$Exp] <- "increase"
set2_fold$Diff[set2_fold$Ctrl > set2_fold$Exp] <- "decrease"
set2_fold$logCtrl <- log2(set2_fold$Ctrl)
set2_fold$logExp <- log2(set2_fold$Exp)
set2_fold$Fold <- set2_fold$logExp - set2_fold$logCtrl




#############
### set 3 ###
#############

# subset set3, day 7
set3_info <- subset(info, Set == "set3" & Day == "day7")
set3_relabun <- subset(relabun, rownames(relabun) %in% set3_info$Sample_ID)
set3_relabun <- set3_relabun[colSums(set3_relabun) > 0]
ncol(set3_relabun) # 1447
set3_relabun[set3_relabun == 0] <- minval


# controls
set3_ctrl_relabun <- set3_relabun[grep("control", rownames(set3_relabun)),]
set3_ctrl.m <- melt(data.frame(SampleID = rownames(set3_ctrl_relabun), set3_ctrl_relabun),
                    value.name = "Ctrl", variable.name = "ASV")
set3_ctrl.m$Rep <- sapply(strsplit(set3_ctrl.m$SampleID, "_"), '[', 4)
set3_ctrl.m$Set <- "set3"
set3_ctrl.m <- set3_ctrl.m[c("Set", "Rep", "ASV", "Ctrl")]


# Na2SiO3
set3_ss_relabun <- set3_relabun[grep("Na2SiO3", rownames(set3_relabun)),]
set3_ss.m <- melt(data.frame(SampleID = rownames(set3_ss_relabun), set3_ss_relabun),
                  value.name = "Exp", variable.name = "ASV")
set3_ss.m$Rep <- sapply(strsplit(set3_ss.m$SampleID, "_"), '[', 4)
set3_ss.m$Inhibitor <- "Na2SiO3"
set3_ss.m$Set <- "set3"
set3_ss.m <- set3_ss.m[c("Set", "Rep", "Inhibitor", "ASV", "Exp")]
identical(set3_ss.m$ASV, set3_ctrl.m$ASV); identical(set3_ss.m$Rep, set3_ctrl.m$Rep)
set3_ss.m$Ctrl <- set3_ctrl.m$Ctrl


# NaPO4
set3_so_relabun <- set3_relabun[grep("NaPO4", rownames(set3_relabun)),]
set3_so.m <- melt(data.frame(SampleID = rownames(set3_so_relabun), set3_so_relabun),
                  value.name = "Exp", variable.name = "ASV")
set3_so.m$Rep <- sapply(strsplit(set3_so.m$SampleID, "_"), '[', 4)
set3_so.m$Inhibitor <- "NaPO4"
set3_so.m$Set <- "set3"
set3_so.m <- set3_so.m[c("Set", "Rep", "Inhibitor", "ASV", "Exp")]
identical(set3_so.m$ASV, set3_ctrl.m$ASV); identical(set3_so.m$Rep, set3_ctrl.m$Rep)
set3_so.m$Ctrl <- set3_ctrl.m$Ctrl


# ZnPO4
set3_zo_relabun <- set3_relabun[grep("ZnPO4", rownames(set3_relabun)),]
set3_zo.m <- melt(data.frame(SampleID = rownames(set3_zo_relabun), set3_zo_relabun),
                  value.name = "Exp", variable.name = "ASV")
set3_zo.m$Rep <- sapply(strsplit(set3_zo.m$SampleID, "_"), '[', 4)
set3_zo.m$Inhibitor <- "ZnPO4"
set3_zo.m$Set <- "set3"
set3_zo.m <- set3_zo.m[c("Set", "Rep", "Inhibitor", "ASV", "Exp")]
identical(set3_zo.m$ASV, set3_ctrl.m$ASV); identical(set3_zo.m$Rep, set3_ctrl.m$Rep)
set3_zo.m$Ctrl <- set3_ctrl.m$Ctrl


# combine
set3_fold <- rbind(set3_ss.m, set3_so.m, set3_zo.m)


# fold change
set3_fold$Diff[set3_fold$Ctrl < set3_fold$Exp] <- "increase"
set3_fold$Diff[set3_fold$Ctrl > set3_fold$Exp] <- "decrease"
set3_fold$logCtrl <- log2(set3_fold$Ctrl)
set3_fold$logExp <- log2(set3_fold$Exp)
set3_fold$Fold <- set3_fold$logExp - set3_fold$logCtrl




###############
### combine ###
###############

folds <- rbind(set1_fold, set2_fold, set3_fold)
folds <- subset(folds, Diff != 0)
#write.csv(folds, "Data/folds_final.csv", na = "", row.names = F)




#################
### corroders ###
#################

chemos <- readxl::read_xlsx("../Data/corroders.xlsx")
taxa$Corroder <- FALSE
taxa$Corroder[taxa$Genus %in% chemos$Genus[is.na(chemos$Genus) == F]] <- TRUE
corroders <- subset(taxa, Corroder == TRUE)


# subset corroders
corroders <- subset(folds, ASV %in% corroders$ASV)
corroders <- merge(corroders, taxa[c("ASV", "Label")], by = "ASV")


# summaries
stats <- 
  cbind(aggregate(Fold ~ Set, mean, data = corroders),
        aggregate(Fold ~ Set, sd, data = corroders)[2])
colnames(stats) <- c("Var", "Mean", "SD")
temp <- 
  cbind(aggregate(Fold ~ Inhibitor, mean, data = corroders),
        aggregate(Fold ~ Inhibitor, sd, data = corroders)[2])
colnames(temp) <- c("Var", "Mean", "SD")
stats <- rbind(stats, temp)
temp <- 
  cbind(aggregate(Fold ~ Label, mean, data = corroders),
        aggregate(Fold ~ Label, sd, data = corroders)[2])
colnames(temp) <- c("Var", "Mean", "SD")
stats <- rbind(stats, temp)
stats$CV <- stats$SD / abs(stats$Mean)


# all
means <- cbind(aggregate(Fold ~ Inhibitor + Set + Label, mean, data = corroders),
               aggregate(Fold ~ Inhibitor + Set + Label, sd, data = corroders)[4])
colnames(means)[4:5] <- c("Mean", "SD")
means$CV <- means$SD / abs(means$Mean)
#write.csv(means, "Data/corroder_stats.csv", row.names = F, na = "")


# plot
inhiblabs <- c(Na2SiO3 = expression(Na[2]*SiO[3]),
               NaPO4 = expression(NaPO[4]),
               ZnPO4 = expression(ZnPO[4]))

means$Label <- gsub(" \\(", "\n(", means$Label)

corr.plot <- 
  ggplot(means, aes(x = Mean, y = Label, fill = Inhibitor, color = Inhibitor)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80") +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), 
           alpha = 0.5, linewidth = 0.3) +
  geom_errorbar(data = means, aes(xmin = Mean - SD, xmax = Mean + SD),
                position = position_dodge2(width = 0.9, preserve = "single", padding = 0.5),
                show.legend = F, size = 0.3) +
  facet_grid(Label ~ Set, scales = "free_y",
             labeller = labeller(Set = c(set1 = "Deep water\nNormal inhibitor levels",
                                         set2 = "Shallow water\nNormal inhibitor levels",
                                         set3 = "Shallow water\nHigh inhibitor levels"))) +
  scale_fill_manual(values = brewer.pal(3, "Dark2"),
                    labels = inhiblabs) +
  scale_color_manual(values = brewer.pal(3, "Dark2"),
                    labels = inhiblabs) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black", face = "italic"),
        axis.title.x = element_text(size = 9, color = "black", face = "bold"),
        axis.title.y = element_text(size = 9, color = "black", face = "bold"),
        legend.title = element_text(size = 9, color = "black", face = "bold"),
        legend.text = element_text(size = 8, hjust = 0),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 9, color = "black", face = "bold"),
        strip.background = element_blank(),
        panel.border = element_rect(linewidth = 0.25, color = "grey80", fill = NA),
        panel.spacing.y = unit(0, "cm"),
        axis.ticks = element_line(linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25)) +
  guides(fill = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
  labs(y = "Phylum (Genus)", 
       x = bquote(bold(Log[2]*" mean fold change from control")))
corr.plot

#ggsave("PLOTS/corroders.pdf", corr.plot, device = "pdf", width = 7.2, height = 5, unit = "in")




#################
### pathogens ###
#################

# subset pathogens
pathogens <- subset(folds, ASV %in% subset(taxa, Genus %in% 
                    c("Pseudomonas", "Legionella", "Mycobacterium", "Arcobacter"))$ASV)
pathogens <- merge(pathogens, taxa[c("ASV", "Label")], by = "ASV")


# summaries
stats <- 
  cbind(aggregate(Fold ~ Set, mean, data = pathogens),
        aggregate(Fold ~ Set, sd, data = pathogens)[2])
colnames(stats) <- c("Var", "Mean", "SD")
temp <- 
  cbind(aggregate(Fold ~ Inhibitor, mean, data = pathogens),
        aggregate(Fold ~ Inhibitor, sd, data = pathogens)[2])
colnames(temp) <- c("Var", "Mean", "SD")
stats <- rbind(stats, temp)
temp <- 
  cbind(aggregate(Fold ~ Label, mean, data = pathogens),
        aggregate(Fold ~ Label, sd, data = pathogens)[2])
colnames(temp) <- c("Var", "Mean", "SD")
stats <- rbind(stats, temp)
stats$CV <- stats$SD / abs(stats$Mean)


# all
means <- cbind(aggregate(Fold ~ Inhibitor + Set + Label, mean, data = pathogens),
               aggregate(Fold ~ Inhibitor + Set + Label, sd, data = pathogens)[4])
colnames(means)[4:5] <- c("Mean", "SD")
means$CV <- means$SD / abs(means$Mean)
#write.csv(means, "Data/pathogen_stats.csv", row.names = F, na = "")



# plot
inhiblabs <- c(Na2SiO3 = expression(Na[2]*SiO[3]),
               NaPO4 = expression(NaPO[4]),
               ZnPO4 = expression(ZnPO[4]))

means$Label <- gsub(" \\(", "\n(", means$Label)

path.plot <- 
  ggplot(means, aes(x = Mean, y = Label, fill = Inhibitor, color = Inhibitor)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80") +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), 
           alpha = 0.5, linewidth = 0.3) +
  geom_errorbar(data = means, aes(xmin = Mean - SD, xmax = Mean + SD),
                position = position_dodge2(width = 0.9, preserve = "single", padding = 0.5),
                show.legend = F, size = 0.3) +
  facet_grid(Label ~ Set, scales = "free_y",
             labeller = labeller(Set = c(set1 = "Deep water\nNormal inhibitor levels",
                                         set2 = "Shallow water\nNormal inhibitor levels",
                                         set3 = "Shallow water\nHigh inhibitor levels"))) +
  scale_fill_manual(values = brewer.pal(3, "Dark2"),
                    labels = inhiblabs) +
  scale_color_manual(values = brewer.pal(3, "Dark2"),
                     labels = inhiblabs) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black", face = "italic"),
        axis.title.x = element_text(size = 9, color = "black", face = "bold"),
        axis.title.y = element_text(size = 9, color = "black", face = "bold"),
        legend.title = element_text(size = 9, color = "black", face = "bold"),
        legend.text = element_text(size = 8, hjust = 0),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 9, color = "black", face = "bold"),
        strip.background = element_blank(),
        panel.border = element_rect(linewidth = 0.25, color = "grey80", fill = NA),
        panel.spacing.y = unit(0, "cm"),
        axis.ticks = element_line(linewidth = 0.25),
        axis.line = element_line(linewidth = 0.25)) +
  guides(fill = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
  labs(y = "Phylum (Genus)", 
       x = bquote(bold(Log[2]*" mean fold change from control")))
path.plot

#ggsave("PLOTS/pathogens.pdf", path.plot, device = "pdf", width = 7, height = 2.3, unit = "in")
