---
title: "Fig.4: Relative Abundances of conjugated BAs in HMP cultures"

# install/load packages
install.packages("data.table")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("extrafont")
pkgs <- c("data.table","tidyverse", "ggplot2", "extrafont")
lapply(pkgs, require, character.only = TRUE)

#load and prepare data table from positive ion LC-MS/MS -- for Fig. 4a and 4b
df.lcms <- fread("06232021_HMPCultures_FCM_quant.csv", check.names = FALSE)
colnames(df.lcms)[1] <- "row.ID"
colnames(df.lcms) <- gsub(" Peak area", "", as.matrix(colnames(df.lcms)))
t_df.lcms <- t(df.lcms)
colnames(t_df.lcms) <- df.lcms$row.ID
t_df.lcms <- t_df.lcms[-(1:3),]
t_df.lcms <- setDT(as.data.frame(t_df.lcms), keep.rownames = "filename")
df.lcms.long <- gather(data = t_df.lcms, key = featureID, value = Abundance, -c(1))
metadata.lcms <- fread("ReDU_metadata_HMPcultures_pos.csv", header=TRUE, sep=",")
data.lcms <- merge(metadata.lcms, df.lcms.long, by.x="filename", by.y="filename")

#filter table to only conjugated BA features and rename feature IDs
plotdata.lcms <- data.lcms %>%
    filter(featureID %in% c("6365","4189","10223","4203","4172","4167","6420","5738","6809","5752","8294","7969","16519","5753","10738","10741","8256","5767","9591","8083","7122","4158","5742","4173","8277","5746","5753","8283","4222"))
plotdata.lcms$featureID <- recode(plotdata.lcms$featureID,"6365" = "Ala-DCA",
                                        "4189" = "Arg-DCA",
                                        "10223" = "Asn-DCA",
                                        "4203" = "Glu-DCA", 
                                        "4172" = "Gln-DCA", 
                                        "4167" = "His-DCA", 
                                        "6420" = "Ile/Leu-DCA",
                                        "5738" = "Lys-DCA", 
                                        "6809" = "Met-DCA",
                                        "5752" = "Phe-DCA",
                                        "8294" = "Ser-DCA",
                                        "7969" = "Thr-DCA",
                                        "16519" = "Trp-DCA",
                                        "5753" = "Tyr-DCA",
                                        "10738" = "Cit-DCA",
                                        "10741" = "Orn-DCA",
                                        "8256" = "Ala-CA",
                                        "5767" = "Arg-CA",
                                        "9591" = "Asn-CA",
                                        "8083" = "Asp-CA",
                                        "7122" = "Cys-CA",
                                        "4158" = "Glu-CA",
                                        "5742" = "Gln-CA",
                                        "4173" = "His-CA",
                                        "8277" = "Ile/Leu-CA",
                                        "5746" = "Lys-CA", 
                                        "5753" = "Phe-CA",
                                        "8283" = "Ser-CA",
                                        "4222" = "Thr-CA")

# add/modify levels and factors
plotdata.lcms$featureID <- gsub("Ile/Leu", "Leu", plotdata.lcms$featureID)
plotdata.lcms$Amino.Acid <- substr(plotdata.lcms$featureID, 1,3)
plotdata.lcms$BAcore <- substr(plotdata.lcms$featureID,5,11)
plotdata.lcms$ATTRIBUTE_sampletype <- recode(plotdata.lcms$ATTRIBUTE_sampletype,"culture_0h" = "t=0h culture",
                                        "culture_72h" = "t=72h culture",
                                        "blank_extraction" = "extraction blank",
                                        "blank_media" = "media blank")
levels(plotdata.lcms$ATTRIBUTE_sampletype) <- gsub(" ", "\n", levels(plotdata.lcms$ATTRIBUTE_sampletype))

# define color scheme
group.colors <- c("Ala"="#C337a9","Arg"="#D88C9A","Asn"="#8E7DBE","Asp"="#8a67ee", "Cit"="#86BBD8", "Cys" = "#2660A4", "Gln"="#002366","Glu"="#3891A6","His"="#84c7bc", "Leu"="#226F54","Lys"="#7FB069", "Met"="#d8f890","Orn"="#BABF95", "Phe"="#FFC53A","Ser"="#F4A261","Thr"="#F9B5AC","Trp"="#EE7674","Tyr"="#8F5D5D")

# generate plot for Fig. 4b
boxplot.4b <- ggplot(plotdata.lcms) +
 aes(x = ATTRIBUTE_sampletype, y = Abundance, colour = as.factor(Amino.Acid)) +
 geom_point(position=position_jitterdodge(), size =0.5) +
 scale_colour_manual(values=group.colors,name = "Amino acid Conjugation") +
 theme_classic() +
 xlab(label = "Sample type") +
 ylab(label = "Peak Area") +
 theme(text = element_text(hjust = 0.5, family = "Myriad Web Pro", size =8), axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size = 8))
boxplot.4b + scale_x_discrete(labels = c("extraction\nblank", "media\nblank", "culture\nt=0h", "culture\nt=72h")) + scale_y_continuous(labels = scales::scientific)
ggsave("Fig4b.pdf")
 
# generate plot for Fig. 4a
plotdata.lcms.noblanks <- plotdata.lcms[!(plotdata.lcms$ATTRIBUTE_Class == "not applicable"),]
boxplot.4a <- ggplot(plotdata.lcms.noblanks) +
 aes(x = ATTRIBUTE_Class, y = Abundance, colour = as.factor(Amino.Acid)) +
 geom_point(position=position_jitterdodge(), size =0.5) +
 scale_colour_manual(values=group.colors,name = "Amino acid Conjugation") +
 theme_classic() +
 xlab(label = "Phylogenetic Class") +
 ylab(label = "Peak Area") +
 theme(text = element_text(hjust = 0.5, family = "Myriad Web Pro", size =8), axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0, size = 8)) +
 facet_wrap(~BAcore, ncol=2)
boxplot.4a + scale_y_continuous(labels = scales::scientific)
ggsave("Fig.4a.pdf")

# load and prepare data table from negative ion LC-IMS-MS -- for Fig. 4c
df.ims <- read.csv("BA-AA HMP 14May2021 Skyline Peak Area Export_noisotopes.csv", stringsAsFactors = FALSE)
metadata.ims <- fread("Metadata_HMP_cultures_IMS_neg.csv", header=TRUE, sep=",")
data.ims <- merge(metadata.ims, df.ims, by.x="Skyline_filename", by.y="File.Name")
data.ims$Molecule <- gsub("Ile", "Leu", data.ims$Molecule)
data.ims$BAcore <- substr(data.ims$Molecule,5,11)
data.ims$Amino.Acid <- substr(data.ims$Molecule, 1,3)

# remove data for unconjugated BAs
data.ims <-data.ims[!(data.ims$Molecule == "CA" | data.ims$Molecule == "CDCA" | data.ims$Molecule == "DCA"| data.ims$Molecule == "HDCA" | data.ims$Molecule == "UDCA"| data.ims$Molecule == "MCA (a)"| data.ims$Molecule == "MCA (b)" | data.ims$Molecule == "MCA (g)"),]

#filter to include only CA, DCA, and CDCA conjugates for plotting (the ones with highest abundance)
plotdata.ims <- data.ims %>%
    filter(BAcore %in% c("CA", "CDCA", "DCA"))

# more filtering for plots -- remove data for L-DOPA conjugate (not observed), remove data with peak areas less than 2000 (LOD for IMS), and remove data for blanks
plotdata.ims <-data[!(plotdata.ims$Amino.Acid == "DOP"),]
plotdata.ims$Area <- as.numeric(plotdata.ims$Area)
plotdata.ims$Area[plotdata.ims$Area < 2000] <- NA
plotdata.ims <- plotdata.ims[complete.cases(plotdata.ims),]
plotdata.ims <- plotdata.ims[!(plotdata.ims$ATTRIBUTE_Genus == "not applicable"),]
write.csv(plotdata.ims, "test.csv")

# order/recode the levels for plotting according to phylogenetic tree
plotdata.ims$ATTRIBUTE_Genus <- factor(plotdata.ims$ATTRIBUTE_Genus, levels = c("Bacteroides","Prevotella","Fusobacterium","Acidaminococcus","Enterococcus", "Lactobacillus", "Clostridium", "unclassified Clostridiales (miscellaneous)", "unclassified Lachnospiraceae", "Cellulosilyticum", "Lachnoclostridium","Catabacter","Peptoniphilus", "Bifidobacterium","Cutibacterium", "Collinsella", "not specified"))
plotdata.ims$ATTRIBUTE_Genus <- recode(plotdata.ims$ATTRIBUTE_Genus, "unclassified Clostridiales (miscellaneous)" = "unclassified Clostridiales")
plotdata.ims <- plotdata.ims %>%
                    filter(!is.na(ATTRIBUTE_Genus))
plotdata.ims$BAcore <- factor(plotdata.ims$BAcore, levels = c("CA","DCA","CDCA"))


# generate plot for Fig. 4c
boxplot.4c <- ggplot(plotdata.ims) +
 aes(x = ATTRIBUTE_Genus, y = as.numeric(Area), colour = Amino.Acid) +
 geom_point(size = 0.65, position=position_jitterdodge(jitter.width=0.1, dodge.width = 1), show.legend = FALSE) +
 theme_classic() +
 scale_colour_manual(values=group.colors,name = "Bile Acid") +
 xlab(label = "Bacterial Genus") +
 ylab(label = "log(Peak Area)") +
 theme(text = element_text(hjust = 0.5, size = 8, family = "Myriad Web Pro"), axis.text.x = element_text(angle = -90, vjust = 0, size = 8)) + 
 facet_wrap(~BAcore, ncol = 4, scales = "free_y")
boxplot.4c + scale_y_continuous(labels = scales::scientific)
ggsave("Fig.4c.pdf")
