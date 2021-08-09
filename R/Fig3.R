---
title: 'Fig3: analysis of iHMP2 metabolomics data'
---
#install/load libraries
# install/load packages
install.packages("data.table", type = "binary")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("extrafont")
install.packages("ggpubr")
install.packages("rstatix")
pkgs <- c("data.table", "tidyverse", "ggplot2", "extrafont", "ggpubr")
lapply(pkgs, require, character.only = TRUE)

# import data
df <- fread("iHMP2data.csv", check.names = FALSE)

# create additional columns
## exclude positive ion data for plotting
df$Amino.Acid <- substr(df$Compound, 1,3)
df$Diagnosis <- factor(df$Diagnosis, levels = c("nonIBD", "CD", "UC"))
df$Diagnosis <- as.character(df$Diagnosis)
plotdata <-df[!(df$Method == "HILIC-pos"),]

# create color scheme
group.fill.colors <- c("nonIBD"="#CADDF3","CD"="#F9B5AC","UC"="#84c7bc")
group.colours <- c("nonIBD"= "#86BBD8", "CD" = "#EE7674", "UC"="#3891A6")

## To create Fig. 3a
# add rows with empty values for better spacing in plot
empties <- data.frame(Compound_orig=unique(plotdata$Compound)[-length(unique(plotdata$Compound))])
empties$Compound <- paste0(empties$Compound_orig,empties$Compound_orig)
empties$value <- NA
empties <- empties %>% add_row(Compound_orig = "Phe-DCA", Compound = "Phe-DCAPhe-DCA", value = "NA")
plotdata_wspaces <- plyr::rbind.fill(plotdata,empties)
plotdata_wspaces$Diagnosis <- factor(plotdata_wspaces$Diagnosis, levels = c("nonIBD", "CD", "UC"))

# generate plot for Fig.3a
boxplot <- ggplot(plotdata_wspaces) +
 aes(x = Compound, y = as.numeric(Abundance), fill = Diagnosis, colour = Diagnosis) +
 geom_boxplot(width = 1, lwd = 0.4, alpha = 0.5, position = position_dodge(width = 1.5), outlier.shape = NA) +
 theme_classic() +
 scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
 scale_colour_manual(values=group.colours,name = "Diagnosis") +
 scale_x_discrete(breaks=unique(df2$Compound)) +
 theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
 xlab(label = "Conjugated Bile Acid") +
 ylab(label = "Peak Area") +
 theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0, size = 8)) 
boxplot + coord_fixed(ratio = 1/200000, ylim = c(0,5000000), expand = FALSE)
ggsave("Fig3a.pdf")

## To create Fig.3b
# separate data by amino acid and BA core
plotdata <- plotdata %>%
  separate(Compound, c("Amino.Acid","BACore"), "-")

# subset data by amino acids and reorder factors
Asn.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Asn"))
Asp.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Asp"))
Cit.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Citrulline"))
Gln.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Gln"))
Glu.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Glu"))
His.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("His"))
Ile.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Ile/Leu"))
Met.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Met"))
Phe.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Phe"))
Thr.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Thr"))
Trp.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Trp"))
Tyr.data <- plotdata %>%
    dplyr::filter(Amino.Acid %in% c("Tyr"))
Asn.data$Diagnosis <- factor(Asn.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Asp.data$Diagnosis <- factor(Asp.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Cit.data$Diagnosis <- factor(Cit.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Gln.data$Diagnosis <- factor(Gln.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Glu.data$Diagnosis <- factor(Glu.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
His.data$Diagnosis <- factor(His.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Ile.data$Diagnosis <- factor(Ile.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Met.data$Diagnosis <- factor(Met.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Phe.data$Diagnosis <- factor(Phe.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Thr.data$Diagnosis <- factor(Thr.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Trp.data$Diagnosis <- factor(Trp.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Tyr.data$Diagnosis <- factor(Tyr.data$Diagnosis, levels = c("nonIBD", "CD", "UC"))
Gln.data$BACore <- factor(Gln.data$BACore, levels = c("CA", "CDCA", "HDCA/UDCA"))
Tyr.data$BACore <- factor(Tyr.data$BACore, levels = c("CA", "CDCA", "HDCA/UDCA"))

# perform pairwise statistical analysis of data using Wilcoxon test and Benjamini-Hochberg correction for p values
stat.test.Asn <- Asn.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Asn <- stat.test.Asn %>% add_xy_position(x = "BACore")

stat.test.Asp <- Asp.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Asp <- stat.test.Asp %>% add_xy_position(x = "BACore")

stat.test.Cit <- Cit.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Cit <- stat.test.Cit %>% add_xy_position(x = "BACore")

stat.test.Gln <- Gln.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Gln <- stat.test.Gln %>% add_xy_position(x = "BACore")

stat.test.Glu <- Glu.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Glu <- stat.test.Glu %>% add_xy_position(x = "BACore")

stat.test.His <- His.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.His <- stat.test.His %>% add_xy_position(x = "BACore")

stat.test.Ile <- Ile.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Ile <- stat.test.Ile %>% add_xy_position(x = "BACore")

stat.test.Met <- Met.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Met <- stat.test.Met %>% add_xy_position(x = "BACore")

stat.test.Phe <- Phe.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Phe <- stat.test.Phe %>% add_xy_position(x = "BACore")

stat.test.Thr <- Thr.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Thr <- stat.test.Thr %>% add_xy_position(x = "BACore")

stat.test.Trp <- Trp.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Trp <- stat.test.Trp %>% add_xy_position(x = "BACore")

stat.test.Tyr <- Tyr.data %>%
  group_by(BACore) %>%
  pairwise_wilcox_test(Abundance_0_adj ~ Diagnosis, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Tyr <- stat.test.Tyr %>% add_xy_position(x = "BACore")

# generate boxplots with adjusted p values
boxplot.Asn.p <- ggboxplot(Asn.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Asn,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Asn.p.pdf", width = 2, height= 2)

boxplot.Asp.p <- ggboxplot(Asp.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Asp,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Asp.p.pdf", width = 2, height= 2)

boxplot.Cit.p <- ggboxplot(Cit.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Cit,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Cit.p.pdf", width = 2, height= 2)

boxplot.Gln.p <- ggboxplot(Gln.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Gln,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Gln.p.pdf", width = 2, height= 2)

boxplot.Glu.p <- ggboxplot(Glu.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Glu,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Glu.p.pdf", width = 2, height= 2)

boxplot.His.p <- ggboxplot(His.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.His,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.His.p.pdf", width = 2, height= 2)

boxplot.Ile.p <- ggboxplot(Ile.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Ile,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Ile.p.pdf", width = 1.1, height= 2)

boxplot.Met.p <- ggboxplot(Met.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Met,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Met.p.pdf", width = 2, height= 2)

boxplot.Phe.p <- ggboxplot(Phe.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Phe,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Phe.p.pdf", width = 2, height= 2)

boxplot.Thr.p <- ggboxplot(Thr.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Thr,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Thr.p.pdf", width = 2, height= 2)

boxplot.Trp.p <- ggboxplot(Trp.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Trp,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Trp.p.pdf", width = 2, height= 2)

boxplot.Tyr.p <- ggboxplot(Tyr.data, "BACore", "Abundance", color = "Diagnosis", fill = "Diagnosis", add = "jitter", add.params = list(size =0.01, jitter = 0.2)) +
theme_classic() +
scale_fill_manual(values=group.fill.colors,name = "Diagnosis") +
stat_pvalue_manual(
    stat.test.Tyr,label = "p.adj", tip.length = 0.01, label.size = 2,
    hide.ns = TRUE
    ) +
scale_colour_manual(values=group.colours,name = "Diagnosis") +
theme(text = element_text(size = 8, family = "Myriad Web Pro")) +
xlab(label = "Conjugated Bile Acid") +
ylab(label = "Peak Area") +
theme(plot.title = element_text(hjust = 0, size = 8), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8)) +
theme(legend.position = "none") +
facet_grid(~Amino.Acid)
ggsave("boxplot.Tyr.p.pdf", width = 2, height= 2)
