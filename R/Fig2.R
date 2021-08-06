---
title: "Fig2: analysis of MASST results"

# install/load packages
install.packages("data.table")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("extrafont")
install.packages("vegan")
pkgs <- c("data.table", "tidyverse", "ggplot2", "extrafont", "vegan")
lapply(pkgs, require, character.only = TRUE)

# combine presence-absence feature table with metadata to obtain data table for analysis
feature.table <- fread('MASST_ConjugatedBAs_FeatureTable.csv', header = TRUE)
metadata <- fread('all_sampleinformation.tsv', header = TRUE)
data <- merge(metadata, feature.table, by.x="sample_name", by.y="filename")

# filter based on inclusion of only data from animals
animal.data <- data %>%
    filter(SampleType %in% c("animal"))

## For Fig.2b -- Tissue/biofluid distribution
# consolidate different types of skin samples to just 'skin'
animal.data <- animal.data %>%
mutate(UBERONBodyPartName = ifelse(UBERONBodyPartName == "skin of body","skin",
                               ifelse(UBERONBodyPartName == "skin of pes", "skin",
                                 ifelse(UBERONBodyPartName == "skin of trunk", "skin",
                                    ifelse(UBERONBodyPartName == "arm skin", "skin",  
                                      ifelse(UBERONBodyPartName == "head or neck skin", "skin",
                                        ifelse(UBERONBodyPartName == "skin of body", "skin",
                                          ifelse(UBERONBodyPartName == "axilla skin", "skin", animal.data$UBERONBodyPartName))))))))

# consolidate different types of blood samples to just 'blood'
animal.data <- animal.data %>%
mutate(UBERONBodyPartName = ifelse(UBERONBodyPartName == "blood serum", "blood",
                                ifelse(UBERONBodyPartName == "blood plasma", "blood", animal.data$UBERONBodyPartName)))

# select metadata category to tabulate -- UBERONBodyPartName
bodypart.data <- aggregate(animal.data[,30:71], list(animal.data$UBERONBodyPartName), sum)
colnames(bodypart.data)[1] <- "UBERON_BodyPart"

# apply function to obtain the proportion of each bile acid detected in different body sites -- each column sums to 1 
bodypart.data[, -1] <- lapply(bodypart.data[ , -1], function(x) x/sum(x, na.rm=TRUE) )

# filter data to include only certain body parts
bodypart.data.filtered <- bodypart.data %>%
    filter(UBERON_BodyPart %in% c("blood", "caecum", "colon", "duodenum", "feces", "gall bladder", "ileum", "jejunum", "liver","oral cavity", "skin", "stomach", "upper digestive tract","urine"))

# transform to long table
## add columns to include number of OH and conjugated amino acid
### change N/As to zeros
bodypart.data.long <- gather(data = bodypart.data.filtered, key = Bile.Acid, value = Count, -c(1))
colnames(bodypart.data.long)[1] <- "UBERON_BodyPart"
bodypart.data.long$Bile.Acid <- gsub("Ile/Leu", "Leu", bodypart.data.long$Bile.Acid)
bodypart.data.long$numberofOH <- substr(bodypart.data.long$Bile.Acid, 7, 7)
bodypart.data.long$aminoacid <- substr(bodypart.data.long$Bile.Acid, 1,3)
bodypart.data.long$aminoacid <- gsub("Leu", "Ile/Leu", bodypart.data.long$aminoacid)
bodypart.data.long <- bodypart.data.long %>%
    mutate(numberofOH = recode(numberofOH,
    "2" = "dihydroxy BA",
    "3" = "trihydroxy BA"))
bodypart.data.long[is.na(bodypart.data.long)] <- 0

# produce figure separated into dihydroxy vs. trihydroxy bile acids -- SI Fig. 20
bodypart.data.heatmap <- ggplot(data = bodypart.data.long, mapping = aes(x = UBERON_BodyPart,
                                                       y = aminoacid,
                                                    fill = as.numeric(Count))) +
    theme_classic() +
    geom_tile(colour = "white") +
    scale_fill_gradient(low = "white", high = "#2660A4", name = "Proportion of matches") +
    xlab(label = "Tissue/biofluid type") +
    ylab(label = "Amino acid conjugation") +
        theme(text = element_text(size = 8, family = "Myriad Web Pro"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), axis.text.y = element_text(size = 8), strip.text.x = element_text(size = 8)) +
    facet_grid(~numberofOH)
bodypart.data.heatmap + coord_fixed()
ggsave("SIFig20.pdf")

# sum values for dihydroxy and trihydroxy bile acids to make combined plot for Fig. 2b
bodypart.data.long.combined <- bodypart.data.long %>%
      group_by(UBERON_BodyPart, aminoacid) %>%
      summarize(Count = sum(Count))
bodypart.data.long.combined['Name']='Tissue/biofluid type'

# generate and plot Fig.2b
bodypart.data.heatmap <- ggplot(data = bodypart.data.long.combined, mapping = aes(x = UBERON_BodyPart,
                                                       y = aminoacid,
                                                    fill = as.numeric(Count))) +
    theme_classic() +
    geom_tile(colour = "white") +
    scale_fill_gradient(low = "white", high = "#2660A4", name = "Proportion of matches") +
    xlab(label = "Tissue/biofluid type") +
    ylab(label = "Amino acid conjugation") +
        theme(text = element_text(size = 8, family = "Myriad Web Pro"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), axis.text.y = element_text(size = 8), strip.text.x = element_text(size = 8)) +
    facet_grid(~Name)
bodypart.data.heatmap + coord_fixed()
ggsave("Fig2b_Bodyparts.pdf")

## For Fig.2b -- Distribution across different species
# select metadata category to tabulate -- NCBITaxonomy
species.data <- aggregate(animal.data[,30:71], list(animal.data$NCBITaxonomy), sum)
colnames(species.data)[1] <- "NCBITaxonomy"

# apply function to obtain the proportion of each bile acid detected in different species -- each column sums to 1 
species.data[, -1] <- lapply(species.data[ , -1], function(x) x/sum(x, na.rm=TRUE) )

# relabel factors to "mouse", "human", "rat", and "other vertebrate"
species.data <- species.data %>%
  mutate(simple.NCBI = ifelse(NCBITaxonomy == "10088|Mus" | NCBITaxonomy == "10090|Mus musculus", "mouse",
                              ifelse(NCBITaxonomy == "10114|Rattus" | NCBITaxonomy == "10116|Rattus norvegicus", "rat",
                                     ifelse(NCBITaxonomy == "9606|Homo sapiens" ,"human",
                                            ifelse(NCBITaxonomy == "not specified" ,"not specified","other vertebrate")))))

#transform to long table
## add columns to include number of OH and conjugated amino acid
### change N/As to zeros
#### order factors for illustration
species.data.long <- gather(data = species.data, key = Bile.Acid, value = Count, -c(1,44))
species.data.long$Bile.Acid <- gsub("Ile/Leu", "Leu", species.data.long$Bile.Acid)
species.data.long$numberofOH <- substr(species.data.long$Bile.Acid, 7, 7)
species.data.long$aminoacid <- substr(species.data.long$Bile.Acid, 1,3)
species.data.long <- species.data.long %>%
  mutate(numberofOH = recode(numberofOH,
                             "2" = "dihydroxy BA",
                             "3" = "trihydroxy BA"))
species.data.long[is.na(species.data.long)] <- 0
species.data.long$simple.NCBI <- factor(species.data.long$simple.NCBI, levels = c("human", "mouse","rat", "other vertebrate", "not specified"))
species.data.long <- species.data.long %>%
  group_by(simple.NCBI, aminoacid) %>%
  summarize(Count = sum(Count))

# generate and plot Fig 2b species distribution heatmap
species.data.long['Name']='Species'
species.data.heatmap <- ggplot(data = species.data.transformed.mod, mapping = aes(x = simple.NCBI,
                                                                                  y = aminoacid,
                                                                                  fill = as.numeric(Count))) +
  theme_classic() +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#2660A4", name = "Proportion of matches") +
  xlab(label = "Species") +
  ylab(label = "Amino acid conjugation") +
  theme(text = element_text(size = 8, family = "Myriad Web Pro"), axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), strip.text.x = element_text(size = 8)) +
  facet_grid(~Name)
species.data.heatmap + coord_fixed()
ggsave("Fig2b_species.pdf")

## For Fig.2c -- Distribution across disease phenotypes
# rename animal samples with DOIDCommonName = "not applicable" to "not specified"
animal.data <- animal.data %>%
mutate(DOIDCommonName = ifelse(DOIDCommonName == "not applicable","not specified", animal.data$DOIDCommonName))

# select metadata category to tabulate -- DOIDCommonName
health.data <- aggregate(animal.data[,30:71], list(animal.data$DOIDCommonName), sum)

# apply function to obtain the proportion of each bile acid detected in different disease phenotypes -- each column sums to 1 
health.data[, -1] <- lapply( health.data[ , -1], function(x) x/sum(x, na.rm=TRUE) )

# filter data to include only those diseases where new conjugated BAs are detected
colnames(health.data)[1] <- "DOIDCommonName"
health.data.filtered <- health.data %>%
  filter(DOIDCommonName %in% c("Chagas disease","circadian rhythm disorders", "Crohn's disease", "disease NOS", "inflammatory bowel disease", "no disease reported", "not specified", "obesity", "sleep deprivation", "ulcerative colitis"))

# transform to long table
## add columns to include number of OH and conjugated amino acid
### change N/As to zeros
health.data.long <- gather(data = health.data.filtered, key = Bile.Acid, value = Count, -c(1))
colnames(health.data.long)[1] <- "DOIDCommonName"
health.data.long$Bile.Acid <- gsub("Ile/Leu", "Leu", health.data.long$Bile.Acid)
health.data.long$numberofOH <- substr(health.data.long$Bile.Acid, 7, 7)
health.data.long$aminoacid <- substr(health.data.long$Bile.Acid, 1,3)
health.data.long <- health.data.long %>%
  mutate(numberofOH = recode(numberofOH,
                             "2" = "dihydroxy BA",
                             "3" = "trihydroxy BA"))
health.data.long[is.na(health.data.long)] <- 0

#recode and reorder levels for illustration
health.data.long$DOIDCommonName <- recode(health.data.long$DOIDCommonName, "Chagas disease" = "Chagas disease",
                                                     "circadian rhythm disorders" = "circadian disorders",
                                                     "sleep deprivation" = "sleep deprivation",
                                                     "Crohn's disease" = "CD", 
                                                     "ulcerative colitis" = "UC", 
                                                     "inflammatory bowel disease" = "IBD", 
                                                     "obesity" = "obesity",
                                                     "disease NOS" = "disease NOS", 
                                                     "not specified" = "not specified", 
                                                     "no disease reported" = "no disease reported")
health.data.long$DOIDCommonName <- factor(health.data.long$DOIDCommonName, levels = c("not specified", "no disease reported","Chagas disease", "circadian disorders", "sleep deprivation","CD", "UC", "IBD", "obesity", "disease NOS"))

health.data.heatmap <- ggplot(data = health.data.transformed.mod, mapping = aes(x = DOIDCommonName,
                                                                                y = aminoacid,
                                                                                fill = as.numeric(Count))) +
  theme_classic() +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "white", high = "#2660A4", name = "Proportion of matches") +
  xlab(label = "Health Condition") +
  ylab(label = "Amino acid conjugation") +
  theme(text = element_text(size = 8, family = "Myriad Web Pro"), axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), strip.text.x = element_text(size = 8)) +
  facet_grid(~numberofOH)
health.data.heatmap + coord_fixed()
ggsave("Fig2c.pdf")

## For Fig.2d -- PCoA using Jaccard distance
# remove Gly and Tau conjugates for analysis
animal.data$Gly_OH3 <- NULL
animal.data$Tau_OH3 <- NULL
animal.data$Gly_OH2 <- NULL
animal.data$Tau_OH2 <- NULL

# remove samples without new conjugated BAs -- aka those that only had Gly or Tau conjugates
plotdata <- animal.data[rowSums(animal.data[, -(1:29)])>0, ]

#remove metadata columns to create matrix
mat = as.matrix(plotdata[,-(1:29)])

# create jaccard distance matrix
jaccard <- vegan::vegdist(mat, method = "jaccard", binary = TRUE, na.rm = FALSE)
jaccard.mat <- as.matrix(jaccard)

# perform PCoA and check goodness of fit 
## determine percentage variance explained for axes
### add metadata back
#### rename columns with PC coordinates
pcoa <- cmdscale(jaccard.mat, k=2, eig = TRUE)
pcoa$GOF
round(pcoa$eig*100/sum(pcoa$eig),1)
pcoa.metadata <- cbind(plotdata,pcoa$points, check.names = TRUE)
names(pcoa.metadata)[68:69] <- c('PC1', 'PC2')

# create groups and colors for visualization
groups <- pcoa.metadata$HealthStatus
pcoa.metadata <- pcoa.metadata %>%
  mutate(HealthStatus = ifelse(HealthStatus == "not applicable","not specified", pcoa.metadata$HealthStatus))
group.fill.colors <- c("chronic illness"="#F36B59","unhealthy (NOS)"="#F9B5AC","healthy"="#54B0A1","not specified" = "#2660A4")

# generate plot -- insert % variance explained for axes
PCoA <- ggplot(pcoa.metadata, aes(x = PC1, y = PC2, colour = HealthStatus)) +
  geom_point(aes(colour = HealthStatus), size =1) +
  scale_colour_manual(name = "Health Status", values=group.fill.colors) +
  theme_classic() +
  xlab(label = "PC1 (15.7%)") +
  ylab(label = "PC2 (12.8%)") +
  theme(text = element_text(size = 8, family = "Myriad Web Pro"))+
  theme(panel.grid.major = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor = element_blank())
PCoA + theme(legend.position = "bottom") + guides(colour=guide_legend(nrow=2, byrow=TRUE))
ggplot("PCoA.pdf")