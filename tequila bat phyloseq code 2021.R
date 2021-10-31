
#       _____                  _  _         ___        _     ___             _           _   		#
#      |_   _|___  __ _  _  _ (_)| | __ _  | _ ) __ _ | |_  | _ \ _ _  ___  (_) ___  __ | |_ 		#
#        | | / -_)/ _` || || || || |/ _` | | _ \/ _` ||  _| |  _/| '_|/ _ \ | |/ -_)/ _||  _|		#
#        |_| \___|\__, | \_,_||_||_|\__,_| |___/\__,_| \__| |_|  |_|  \___/_/ |\___|\__| \__|		#
#                    |_|                                                  |__/              		  #



#Created by Luis Víquez-R, University of Ulm, Germany
#contact me at luis.viquez@alumni.uni-ulm.de

                                          
#install.packages("vegan")
#install.packages("tidyverse")
#install.packages("ggplot2")

#R version 3.6.0
#install.packages("multcompView")
#install.packages("stringr")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("phyloseq")
#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")
library(metagMisc)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(multcompView)
library(stringr)
library(vegan)
library(btools)

#setwd("D:/migration138a/phyloseq May 2021")
setwd("C:/Users/Luis/Google Drive/PHD/A Faithful Gut/phyloseq May 2021")

#####import phyloseq########
lepto_map <-import_qiime_sample_data ("phyloseq_req/Mapfile20210519.txt")
lepto_tree <- read_tree_greengenes ("phyloseq_req/tree138-2105.nwk")
lepto_biom <-import_biom("phyloseq_req/leptoTaxonomy138-0521.biom")
lepto_biom
head(otu_table(lepto_biom))
head(tax_table(lepto_biom))
head(lepto_map)
#merge into new object
migration<- merge_phyloseq (lepto_biom, lepto_map, lepto_tree)
migration
sample_data(migration)

#quickly check the object structure
otu_table(migration)
tax_table(migration)[1:20,1:7]
phy_tree(migration)
head(sample_data(migration))

####### to change the name of the rows in the tax_table ###### 
yearchange<- data.frame(sample_data(migration, errorIfNULL = FALSE))

yearchange$year<-as.character(yearchange$year)
yearchange[is.na(yearchange)] <- ""
sample_data(migration) <- as.data.frame(yearchange)

#remove un-used IDs
migration <- migration %>%
  subset_samples(Biome!="remove")

# Filter ASVs present in controls from entire dataset
controls <- migration %>%
  subset_samples(species == "control") 
controls@sam_data
controls_1 <- prune_taxa(taxa_sums(controls) > 2, controls)
badtaxa<-taxa_names(controls_1)
alltaxa<-taxa_names(migration)
alltaxa1 <- alltaxa[!(alltaxa %in% badtaxa)]

migration = prune_taxa(alltaxa1, migration)

#Change label for taxonomy
tax <- data.frame(tax_table(migration))
tax
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "d__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""
tax.clean[1:20,1:7]
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

####### Filling the NA's in the tax table#######
tax.clean[is.na(tax.clean)] <- ""

#### re import as matrix into the S4 object
tax_table(migration) <- as.matrix(tax.clean)

tax_table(migration)[1:100,1:7]

#change "NA" for ""
migration = subset_taxa(migration, Phylum!="")
migration

#Dataframe for ASVs counts
sample_sum_df0<- data.frame(sum = sample_sums(migration))
sample_sum_df0

#plot for ASV counts
asvcounts0<-ggplot(sample_sum_df0, aes(x = sum)) + geom_histogram()
asvcounts0+
  ggtitle("Distribution of sample sequencing depth")+
  xlab("Read counts")+
  ylab("Frequency")+theme_classic(base_size = 16)
rank_names(migration)
#remove samples not included in this dataset
lepto_sub <- migration %>%
  subset_samples(Biome!="other")
lepto_sub
head(sample_data(lepto_sub))
migration

#remove samples with less than 8000 reads
lepto_sub<-prune_samples(sample_sums(lepto_sub)>=8000,lepto_sub)
#remove outlier with more than 100000 reads
lepto_sub<-prune_samples(sample_sums(lepto_sub)<=100000,lepto_sub)
sample_sum_df<-data.frame(sum = sample_sums(lepto_sub))
summary(sample_sum_df)
lepto_sub

#basic stats
std_coverage<-sd(sample_sum_df$sum)
mean_coverage<-mean(sample_sum_df$sum)
min_coverage<-min(sample_sum_df$sum)
max_coverage<-max(sample_sum_df$sum)
std_coverage
mean_coverage
min_coverage
max_coverage
#histogram
asvcounts<-ggplot(sample_sum_df, aes(x = sum)) + geom_histogram()
asvcounts+
  ggtitle("Distribution of sample sequencing depth")+
  xlab("Read counts")+
  ylab("Frequency")+theme_classic(base_size = 16)
rank_names(lepto_sub)

################################################################################################################
#####Calculating Alpha diversity Measurements#######
lepto_subr<-rarefy_even_depth(lepto_sub, sample.size = min(sample_sums(lepto_sub)),
                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
min(sample_sums(lepto_sub))
Richness<-estimate_richness(lepto_subr,measures=c("Observed","Shannon"))
FaithsPD<-estimate_pd(lepto_subr)
names(Richness)
sample_data(lepto_subr)
Alpha<-sample_data(lepto_subr)
names(Alpha)

Alpha$Observed<-Richness$Observed
Alpha$Shannon<-Richness$Shannon
Alpha$FPD<-FaithsPD$PD
head(Alpha)

#Lets also add sequencing depth per sample:
Alpha$SequencingDepth<-sample_sums(lepto_subr)
head(Alpha)

###################################ALpha Diversity#####################################################################

###########color pallets##########
print("Color palletes were generated using the IwantHUE online tool, check it out at http://medialab.github.io/iwanthue/")
colas<-c("Chamela"="#9e71c0",
  "Coquimatlán"="#92a619",
  "Pinacate"="#ff5c27",
  "Desert"="#ff5c27",
  "Dry Forest"="#8B7E62")
colastemp<-c("Chamela 2016"="#9e71c0",
  "Chamela 2017"="#6372ad",
  "Pinacate 2015"="#d6a336",
  "Pinacate 2016"="#ef0033",
  "Pinacate 2017"="#ff5c27",
  "Coquimatlán 2017"="#92a619")
colastempbeta<-c("Chamela 2016"="#9e71c0",
             "Chamela 2017"="#6372ad",
             "Pinacate 2015"="#d6a336",
             "Pinacate 2016"="#ef0033",
             "Pinacate 2017"="#ff5c27")
colasbeta<-c("Chamela"="#9e71c0",
         "Coquimatlán"="#92a619",
         "Pinacate"="#ff5c27")
Taxacolors<-c("Streptococcaceae"="#d8ca48",
             "Helicobacteraceae"="#ff759a",
             "Erysipelotrichaceae"="#8edd49",
             "Aeromonadaceae"="#9172e9",
             "Yersiniaceae"="#e54e96",
             "Leuconostocaceae"="#749d40",
             "Lachnospiraceae"="#7689d5",
             "Mycobacteriaceae"="#e75b32",
             "Mycoplasmataceae"="#c449d1",
             "Enterobacteriaceae"="#50adeb",
             "Clostridiaceae"="#a38d2e",
             "Bacillaceae"="#7c79d0",
             "Enterococcaceae"="#b3255c",
             "Actinobacteriota"="#5daac4",
             "Proteobacteria"="#cd8f48",
             "Firmicutes"="#aca5c9",
             "Fusobacteriota"="#7c79d0",
             "Verrucomicrobiota"="#d8ca48",
             "Campilobacterota"="#ceda9e",
             "Methylophilaceae"="#d082be",
             "below cutoff"="grey",
             "Others (< 5% abund.)"="#898181")
years<-c("2015"="#d6a336",
  "2016"="#ef0033",
  "2017"="#ff7851")
ANCOMcolors=c("Bacillaceae"="#4faf5d",
              "Clostridiaceae"="#b648a2",
              "Erysipelotrichaceae"="#88ac3a",
              "Geodermatophilaceae"="#6e5ec2",
              "Helicobacteraceae"="#c2993c",
              "Lachnospiraceae"="#6597bd",
              "Leuconostocaceae"="#c05831",
              "Mycobacteriaceae"="#73aa90",
              "Peptostreptococcaceae"="#bd435c",
              "Rhizobiaceae"="#545b34",
              "Staphylococcaceae"="#c68fbb",
              "Streptococcaceae"="#bc8d78",
              "Corynebacteriaceae"="#aca5c9",
              "Mycoplasmataceae"="#45998e",
              "below cutoff"="grey")
Taxacolors22<-c("Streptococcaceae"="#d8ca48",
              "Helicobacteraceae"="#ff759a",
              "Aeromonadaceae"="#9172e9",
              "Leuconostocaceae"="#749d40",
              "Mycoplasmataceae"="#c449d1",
              "Enterobacteriaceae"="#50adeb",
              "Clostridiaceae"="#a38d2e",
              "Enterococcaceae"="#b3255c",
              "Others (< 5% abund.)"="#898181")

Taxaphyla22<-c("Proteobacteria"="#cd8f48",
         "Firmicutes"="#aca5c9",
         "Actinobacteriota"="#5daac4",
         "Campilobacterota"="#ceda9e",
         "Others (< 5% abund.)"="#898181")

#####plot themes######

themeboxplots<-theme_minimal()+theme(
  axis.title.x = element_text(color = "black", face="bold", size = 16),
  axis.title.y = element_text(color="black",face="bold", size=16),  
  axis.text.x=element_text(colour="black", face="bold", size=14,angle = 0, hjust = 0.5),
  axis.text.y=element_text(colour="black", size = 14),
  title = element_text(color = "black", face="bold", size = 16),
  legend.position='none',
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank())
  
themecomposition<-theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=14),  
        axis.text.x=element_text(colour="black", face="bold", size=14, hjust = 0.5),
        axis.text.y=element_text(colour="black", size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position='right',
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

themebetaplots<-theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.title.y = element_text(face="bold", size=14),  
        axis.text.x=element_text(colour="black", face="bold", size=14, hjust = 0.5),
        axis.text.y=element_text(colour="black", size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position='top',
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

themevolcano<-theme_minimal() +
  theme(text=element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=14),
        axis.title.y = element_text(face="bold", size=14),  
        axis.text.x=element_text(colour="black", face="bold", size=14, hjust = 0.5),
        axis.text.y=element_text(colour="black", size = 14),
        strip.text.x=element_text(face="bold", size=14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position='right',
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.background  = element_rect())


#####All Localites######

#Shannon all sites
Shanbloc<-ggplot(data=Alpha,aes(x=Locality, y=Shannon, fill=Locality))+ggtitle("Shannon diversity per locality")
Shanbloc+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Shannon Index")+scale_fill_manual(values=colas)+themeboxplots
ggsave("news/Shannon boxplot all sites.png",dpi = 300, units = "in", height = 5, width = 4)

#Observed ASVs all sites
ASVbloc<-ggplot(data=Alpha,aes(x=Locality, y=Observed, fill=Locality))+scale_fill_manual(values=colas)+ggtitle("Observed ASVs per locality")
ASVbloc+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Number of ASVs")+themeboxplots
ggsave("news/Observed ASVs in all sites.png",dpi = 300, units = "in", height = 5, width = 4)

#Faith's PD all sites
FPDbloc<-ggplot(data=Alpha,aes(x=Locality, y=FPD, fill=Locality))+ggtitle("Faith's PD per locality")
FPDbloc+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Faith's PD")+scale_fill_manual(values=colas)+themeboxplots
ggsave("news/Faiths PD boxplot all sites.png",dpi = 300, units = "in", height = 5, width = 4)

#### BIOME ######

#Shannon all sites
BiomeShanbloc<-ggplot(data=Alpha,aes(x=Biome, y=Shannon, fill=Biome))+ggtitle("Shannon diversity per Biome")
BiomeShanbloc+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Shannon Index")+scale_fill_manual(values=colas)+themeboxplots
ggsave("news/BiomeShannon boxplot all sites.png",dpi = 300, units = "in", height = 5, width = 4)

#Observed ASVs all sites
BiomeASVbloc<-ggplot(data=Alpha,aes(x=Biome, y=Observed, fill=Biome))+scale_fill_manual(values=colas)+ggtitle("Observed ASVs per Biome")
BiomeASVbloc+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Number of ASVs")+themeboxplots
ggsave("news/BiomeObserved ASVs in all sites.png",dpi = 300, units = "in", height = 5, width = 4)

#Faith's PD all sites
FPDbloc<-ggplot(data=Alpha,aes(x=Biome, y=FPD, fill=Biome))+ggtitle("Faith's PD per Biome")
FPDbloc+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Faith's PD")+scale_fill_manual(values=colas)+themeboxplots
ggsave("news/BiomeFaiths PD boxplot all sites.png",dpi = 300, units = "in", height = 5, width = 4)

###### For Betacross#####
#Subsetting the database using Crossover column that groups Pinacate (2015-2016-2017) and Chamela (2016-2017)
Betacross<- Alpha %>%
  subset_samples(Crossover=="TRUE")

#Shannon Pinacate-Chamela
BetacrossPinChashannon<-ggplot(data=Betacross,aes(x=Betaplus, y=Shannon, fill=Betaplus))+scale_fill_manual(values=colastemp)+ggtitle("Shannon index by Betaplus")
BetacrossPinChashannon+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Shannon Index")+themeboxplots
ggsave("news/Females by betaplus SHANNON.png",dpi = 300, units = "in", height = 4, width = 12)

#ASVs Pinacate-Chamela
BetacrossPinChaASV<-ggplot(data=Betacross,aes(x=Betaplus, y=Observed, fill=Betaplus))+scale_fill_manual(values=colastemp)+ggtitle("ASVs by Betaplus")
BetacrossPinChaASV+geom_boxplot()+geom_jitter(width = 0.25)+ylab("number of ASVs")+themeboxplots
ggsave("news/Females by betaplus ASV.png",dpi = 300, units = "in", height = 4, width = 12)

#DPD Pinacate-Chamela
BetacrossPinChaFPD<-ggplot(data=Betacross,aes(x=Betaplus, y=FPD, fill=Betaplus))+scale_fill_manual(values=colastemp)+ggtitle("ASVs by Betaplus")
BetacrossPinChaFPD+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Faith's PD")+themeboxplots
ggsave("news/Females by betaplus FPD.png",dpi = 300, units = "in", height = 4, width = 12)

####Composition by Locality######

melt_Family <- lepto_sub %>%
  tax_glom(taxrank = "Family") %>%      
  psmelt() %>%                         
  arrange(Family)                     

melt_Family2 <- aggregate(Abundance ~ Family + Locality, 
                          data= melt_Family, 
                          sum)

melt_Family2$Family <- as.character(melt_Family2$Family) 

melt_Family3 <- melt_Family2 %>% 
  dplyr::group_by(Locality) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

#Get Remainers to add in graph
remainers_Family <- melt_Family2 %>% 
  dplyr::group_by(Locality) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Family2 <- aggregate(rel.freq ~ Locality, 
                               data= remainers_Family, 
                               sum)

remainers_Family2$Family <- "Others (< 5% abund.)"

join_Family <- full_join(melt_Family3,remainers_Family2)

join_Family <- join_Family %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Family$Family <- as.factor(join_Family$Family)
join_Family$Family <- reorder(join_Family$Family, join_Family$Abundance)
family_species <- ggplot(join_Family, aes(x = Locality, y = rel.freq, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") + scale_fill_manual(values = Taxacolors)+
  xlab("Locality") +ggtitle("composition by Locality")+
  geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)+themecomposition 
family_species
ggsave("news/Overall Locality composition FAMILY.png",dpi = 300, units = "in", height = 5, width = 10)

####Composition by Biome######

melt_Family <- lepto_sub %>%
  tax_glom(taxrank = "Family") %>%      
  psmelt() %>%                         
  arrange(Family)                     

melt_Family2 <- aggregate(Abundance ~ Family + Biome, 
                          data= melt_Family, 
                          sum)

melt_Family2$Family <- as.character(melt_Family2$Family) 

melt_Family3 <- melt_Family2 %>% 
  dplyr::group_by(Biome) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

#Get Remainers to add in graph
remainers_Family <- melt_Family2 %>% 
  dplyr::group_by(Biome) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Family2 <- aggregate(rel.freq ~ Biome, 
                               data= remainers_Family, 
                               sum)

remainers_Family2$Family <- "Others (< 5% abund.)"

join_Family <- full_join(melt_Family3,remainers_Family2)

join_Family <- join_Family %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Family$Family <- as.factor(join_Family$Family)
join_Family$Family <- reorder(join_Family$Family, join_Family$Abundance)
family_species <- ggplot(join_Family, aes(x = Biome, y = rel.freq, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") + scale_fill_manual(values = Taxacolors)+
  xlab("Biome") +ggtitle("composition by Biome")+
  geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)+
  themecomposition
family_species
ggsave("news/Overall Biome composition FAMILY.png",dpi = 300, units = "in", height = 5, width = 10)

#### composition in Pinacate 2015-2016-2017###########
pinalepto<-lepto_sub %>%
  subset_samples(Crossover==TRUE) 
melt_Familypina <- pinalepto %>%
  tax_glom(taxrank = "Family") %>%      
  psmelt() %>%                           
  arrange(Family)                     

melt_Family2pina <- aggregate(Abundance ~ Family + Betaplus, 
                              data= melt_Familypina, 
                              sum)

melt_Family2pina$Family <- as.character(melt_Family2pina$Family) 

melt_Family3pina <- melt_Family2pina %>% 
  dplyr::group_by(Betaplus) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

#Get Remainers to add in graph
remainers_Familypina <- melt_Family2pina %>% 
  dplyr::group_by(Betaplus) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Family2pina <- aggregate(rel.freq ~ Betaplus, 
                                   data= remainers_Familypina, 
                                   sum)

remainers_Family2pina$Family <- "Others (< 5% abund.)"

join_Familypina <- full_join(melt_Family3pina,remainers_Family2pina)

join_Familypina <- join_Familypina %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Familypina$Family <- as.factor(join_Familypina$Family)
join_Familypina$Family <- reorder(join_Familypina$Family, join_Familypina$Abundance)
family_speciespina<- ggplot(join_Familypina, aes(x = Betaplus, y = rel.freq, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") + ggtitle("Composition by year and locality")+
  ylab("Relative abundance") + scale_fill_manual(values = Taxacolors22)+
  xlab("Year") +
  geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)+
  themecomposition
family_speciespina
ggsave("news/Betaplus composition FAMILY.png",dpi = 600, units = c("in"), height = 6, width = 14)

####Composition by Locality phylum######

melt_Phylum <- lepto_sub %>%
  tax_glom(taxrank = "Phylum") %>%      
  psmelt() %>%                         
  arrange(Phylum)                     

melt_Phylum2 <- aggregate(Abundance ~ Phylum + Locality, 
                          data= melt_Phylum, 
                          sum)

melt_Phylum2$Phylum <- as.character(melt_Phylum2$Phylum) 

melt_Phylum3 <- melt_Phylum2 %>% 
  dplyr::group_by(Locality) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

#Get Remainers to add in graph
remainers_Phylum <- melt_Phylum2 %>% 
  dplyr::group_by(Locality) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Phylum2 <- aggregate(rel.freq ~ Locality, 
                               data= remainers_Phylum, 
                               sum)

remainers_Phylum2$Phylum <- "Others (< 5% abund.)"

join_Phylum <- full_join(melt_Phylum3,remainers_Phylum2)

join_Phylum <- join_Phylum %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Phylum$Phylum <- as.factor(join_Phylum$Phylum)
join_Phylum$Phylum <- reorder(join_Phylum$Phylum, join_Phylum$Abundance)
family_species <- ggplot(join_Phylum, aes(x = Locality, y = rel.freq, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") + scale_fill_manual(values = Taxaphyla22)+
  xlab("Locality") +ggtitle("composition by Locality")+
  geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)+
  themecomposition
family_species
ggsave("news/Overall Locality composition PHYLUM.png",dpi = 600, units = "in", height = 5, width = 14)

#### Betaplus COMPOSITION BY YEARS#####
pinalepto<-lepto_sub %>%
  subset_samples(Crossover==TRUE) 

pinalepto
melt_Phylumpina <- pinalepto %>%
  tax_glom(taxrank = "Phylum") %>%      
  psmelt() %>%                           
  arrange(Phylum)                     

melt_Phylum2pina <- aggregate(Abundance ~ Phylum + Betaplus, 
                              data= melt_Phylumpina, 
                              sum)

melt_Phylum2pina$Phylum <- as.character(melt_Phylum2pina$Phylum) 

melt_Phylum3pina <- melt_Phylum2pina %>% 
  dplyr::group_by(Betaplus) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq > 0.05)

#Get Remainers to add in graph
remainers_Phylumpina <- melt_Phylum2pina %>% 
  dplyr::group_by(Betaplus) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),2)) %>%
  dplyr::filter(rel.freq < 0.05)

remainers_Phylum2pina <- aggregate(rel.freq ~ Betaplus, 
                                   data= remainers_Phylumpina, 
                                   sum)

remainers_Phylum2pina$Phylum <- "Others (< 5% abund.)"

join_Phylumpina <- full_join(melt_Phylum3pina,remainers_Phylum2pina)

join_Phylumpina <- join_Phylumpina %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Phylumpina$Phylum <- as.factor(join_Phylumpina$Phylum)
join_Phylumpina$Phylum <- reorder(join_Phylumpina$Phylum, join_Phylumpina$Abundance)
phylum_speciespina<- ggplot(join_Phylumpina, aes(x = Betaplus, y = rel.freq, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") + ggtitle("Composition by year and locality")+
  ylab("Relative abundance") + scale_fill_manual(values = Taxaphyla22)+
  xlab("Year") +
  geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=7)+
  themecomposition
phylum_speciespina
ggsave("news/Beta plus composition PHYLUM.png",dpi = 600, units = c("in"), height = 6, width = 14)

#generalized linear models for alpha diversity
Alphadf<-data.frame(sample_data(Alpha))
Alphadf$logshannon<-log(Alphadf$Shannon)
Alphadf$logobserved<-(Alphadf$Observed)
Alphadf$sqrtshannon<-sqrt(Alphadf$Shannon)
Alphadf$logFPD<-log(Alphadf$FPD)

#lm locality for FPD
localitylogFPDs<-lm(logFPD~Locality+sex+year, data=Alphadf)
summary(localitylogFPDs)
anova(localitylogFPDs)
#lm for locality for ASV
localitylogObserved<-lm(logobserved~Locality+sex+year, data=Alphadf)
summary(localitylogObserved)
anova(localitylogObserved)
#lm for locality for Shannon
localityShannonsqrt<-lm(sqrtshannon~Locality+sex+year, data=Alphadf)
summary(localityShannonsqrt)
anova(localityShannonsqrt)

######ALPHA DIVERSITY MODELS FOR PINACATE- CHAMELA#####
alphaCross<-Alphadf%>%
  subset(Crossover=="TRUE")
head(alphaCross)

logFPDsCross<-lm(logFPD~year+Biome, data=alphaCross)
summary(logFPDsCross)
anova(logFPDsCross)
LFP<-aov(logFPDsCross)
TukeyHSD(LFP)

logObservedCross<-lm(logobserved~year+Biome, data=alphaCross)
summary(logObservedCross)
anova(logObservedCross)
LOP<-aov(logObservedCross)
TukeyHSD(LOP)

ShannonsqrtCross<-lm(sqrtshannon~year+Biome, data=alphaCross)
summary(ShannonsqrtCross)
anova(ShannonsqrtCross)
SSP<-aov(ShannonsqrtCross)
TukeyHSD(SSP)


#######################BETA DIVERSITY#########################################
###preproccess Beta Diversity####
#remove ALL singletons
lepto_sub1.1<- filter_taxa(lepto_subr, function (x) {sum(x > 0) >1}, prune=TRUE)
lepto_sub1.1
otu_table(lepto_sub1.1)
prevalencedf = apply(X = otu_table(lepto_sub1.1),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

prevalencedf = data.frame(Prevalence = prevalencedf,
                          TotalAbundance = taxa_sums(lepto_sub1.1),
                          tax_table(lepto_sub1.1))
head(prevalencedf)

#we define a a prevalence threshold of 7% of the samples

nsamples(lepto_sub1.1)
prevalenceThreshold = 0.07* nsamples(lepto_sub1.1)
prevalenceThreshold
lepto_sub1.1<- subset_taxa(lepto_sub1.1, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
lepto_sub1.1
plyr::ddply(prevalencedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

prevalencedf[1:10,]
prevalencedf1 = subset(prevalencedf, Phylum %in% get_taxa_unique(lepto_sub1.1, taxonomic.rank = "Phylum"))
prevalencedf1[1:10,]
ggplot(prevalencedf1, aes(TotalAbundance, Prevalence / nsamples(lepto_sub1.1),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
str(prevalencedf)
str(prevalencedf1)

#Now we use this value to filter out any otu that is not present in at least 10% of all samples
#(prevalencedf1$Prevalence >= prevalenceThreshold)
keepTaxa = rownames(prevalencedf1)[(prevalencedf1$Prevalence >= prevalenceThreshold)]
length(keepTaxa)
lepto_sub1 = prune_taxa(keepTaxa, lepto_sub1.1)
lepto_sub1
# weighted unifrac 
DistW = distance(lepto_sub1,method="wunifrac")
#unweighted unifrac 
DistUW = distance(lepto_sub1,method="uunifrac")
#create ordination
ordW = ordinate(lepto_sub1, method = "PCoA", distance = DistW)
ordUW = ordinate(lepto_sub1, method = "PCoA", distance = DistUW)
#lets visualise how informative each ordination is for each distance matrix
plot_scree(ordW)
plot_scree(ordUW)
#plot ordinations
BetaWEallsites<-plot_ordination(lepto_sub1, ordW, color = "Locality",axes=c(2,1))+  stat_ellipse() +  
  ggtitle("Weighted Unifrac by locality")+scale_color_manual(values=colas)+themebetaplots
BetaWEallsites
ggsave("news/PCoA Weighted UniFrac by locality.png",dpi = 300, units = c("in"), height = 8, width = 10)
BetaUWallsites<-plot_ordination(lepto_sub1, ordUW, color = "Locality",axes=c(2,1))+  stat_ellipse() +  
  ggtitle("UnWeighted UniFrac by locality")+scale_color_manual(values=colas)+themebetaplots
BetaUWallsites
ggsave("news/PCoA Un-Weighted UniFrac by locality.png",dpi = 300, units = c("in"), height = 8, width = 10)
#Centroid plots
#### centroids plot for UW Unifrac localities#####
metadata<-data.frame(sample_data(lepto_sub1))
UW_ab<-data.frame(ordUW$vectors[,1],
                  ordUW$vectors[,2])
colnames(UW_ab)[1]<-"MDS1"
colnames(UW_ab)[2]<-"MDS2"
UWcentroid<-cbind(UW_ab,metadata)
centroids_UW <- as.data.frame(UWcentroid %>% 
                                dplyr::group_by(Locality) %>% # calculate functions below for each group
                                dplyr::summarise(mean_MDS1=mean(MDS1),
                                                 mean_MDS2=mean(MDS2),
                                                 n_MDS1=length(MDS1),
                                                 n_MDS2=length(MDS2),
                                                 stdv_MDS1=sd(MDS1),
                                                 stdv_MDS2=sd(MDS2),
                                                 se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                 se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
UW_unifrac_plot_by_Locality <- ggplot(data = centroids_UW, aes(x=mean_MDS1, y=mean_MDS2, color=Locality))+
  geom_point(size=3) +
  geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                    ymax=mean_MDS2+se_MDS2), width=0.02,  size=0.3)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.02, size=0.3) +
  geom_point(data = UWcentroid, aes(x=MDS1, y=MDS2, color=Locality), alpha=0.2, size=4) +
  scale_color_manual(values = colasbeta)+
  labs(x = "MDS1 [16.2%]", y="MDS2 [12.4%]")+
  ggtitle("Unweighted UniFrac") +
  themebetaplots
UW_unifrac_plot_by_Locality+stat_ellipse(data=UWcentroid, aes(x=MDS1, y=MDS2, color=Locality),inherit.aes = FALSE)
ggsave("news/UW_unifrac_plot_by_Locality centr.png",dpi = 300, units = c("in"), height = 8, width = 10)

### Centroid plot UW for localities#####
W_ab<-data.frame(ordW$vectors[,1],
                 ordW$vectors[,2])
colnames(W_ab)[1]<-"MDS1"
colnames(W_ab)[2]<-"MDS2"
Wcentroid<-cbind(W_ab,metadata)
centroids_W <- as.data.frame(Wcentroid %>% 
                               dplyr::group_by(Locality) %>% # calculate functions below for each group
                               dplyr::summarise(mean_MDS1=mean(MDS1),
                                                mean_MDS2=mean(MDS2),
                                                n_MDS1=length(MDS1),
                                                n_MDS2=length(MDS2),
                                                stdv_MDS1=sd(MDS1),
                                                stdv_MDS2=sd(MDS2),
                                                se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
W_unifrac_plot_by_Locality <- ggplot(data = centroids_W, aes(x=mean_MDS1, y=mean_MDS2, color=Locality))+
  geom_point(size=3) +
  geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                    ymax=mean_MDS2+se_MDS2), width=0.03,  size=0.5)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.03, size=0.5) +
  geom_point(data = Wcentroid, aes(x=MDS1, y=MDS2, color=Locality), alpha=0.2, size=4) +
  scale_color_manual(values = colasbeta)+
  labs(x = "MDS1 [28.3%]", y="MDS2 [20.5%]")+
  ggtitle("Weighted UniFrac") +
  themebetaplots
W_unifrac_plot_by_Locality+stat_ellipse(data=Wcentroid, aes(x=MDS1, y=MDS2, color=Locality),inherit.aes = FALSE)
ggsave("news/W_unifrac_plot_by_Locality centr.png",dpi = 300, units = c("in"), height = 8, width = 10)


## BETA JUST Pinacate-Chamela
lepto_sub_PIN_CHA<- lepto_subr %>%
  subset_samples(Crossover=="TRUE")
lepto_sub_PIN_CHA<- lepto_sub_PIN_CHA %>%
  subset_samples(Betaplus!="Coquimatlán 2017")
head(sample_data(lepto_sub_PIN_CHA))
prevalencedf_PIN_CHA = apply(X = otu_table(lepto_sub_PIN_CHA),
                             MARGIN = 1,
                             FUN = function(x){sum(x > 0)})
prevalencedf_PIN_CHA = data.frame(Prevalence = prevalencedf_PIN_CHA,
                                  TotalAbundance = taxa_sums(lepto_sub_PIN_CHA),
                                  tax_table(lepto_sub_PIN_CHA))
head(prevalencedf_PIN_CHA)

#we define a a prevalence threshold of 7% of the samples
prevalenceThreshold_PIN_CHA = 0.07 * nsamples(lepto_sub_PIN_CHA)
prevalenceThreshold_PIN_CHA
lepto_sub_PIN_CHA<- subset_taxa(lepto_sub_PIN_CHA, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
lepto_sub_PIN_CHA
plyr::ddply(prevalencedf_PIN_CHA, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})
prevalencedf_PIN_CHA[1:10,]
prevalencedf1_PIN_CHA = subset(prevalencedf_PIN_CHA, Phylum %in% get_taxa_unique(lepto_sub1.1, taxonomic.rank = "Phylum"))
prevalencedf1_PIN_CHA[1:10,]
ggplot(prevalencedf1_PIN_CHA, aes(TotalAbundance, Prevalence / nsamples(lepto_sub_PIN_CHA),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.5) +
  scale_x_log10() +  xlab("Total Abundance in pinacate") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
str(prevalencedf_PIN_CHA)
str(prevalencedf1_PIN_CHA)

#Now we use this value to filter out any otu that is not present in at least 30% of all samples
(prevalencedf1_PIN_CHA$Prevalence >= prevalenceThreshold_PIN_CHA)
keepTaxa_PIN_CHA = rownames(prevalencedf1_PIN_CHA)[(prevalencedf1_PIN_CHA$Prevalence >= prevalenceThreshold_PIN_CHA)]
length(keepTaxa_PIN_CHA)
lepto_sub1_PIN_CHA = prune_taxa(keepTaxa_PIN_CHA, lepto_sub_PIN_CHA)
lepto_sub1_PIN_CHA

#weighted unifrac
DistW_PIN_CHA= distance(lepto_sub1_PIN_CHA,method="wunifrac")
#Unweighted unifrac
DistUW_PIN_CHA= distance(lepto_sub1_PIN_CHA,method="uunifrac")
#create ordination
ordW_PIN_CHA = ordinate(lepto_sub1_PIN_CHA, method = "PCoA", distance = DistW_PIN_CHA)
ordUW_PIN_CHA= ordinate(lepto_sub1_PIN_CHA, method = "PCoA", distance = DistUW_PIN_CHA)
BetaWpinacate_Chamela<-plot_ordination(lepto_sub1_PIN_CHA, ordW_PIN_CHA, color = "Betaplus",axes=c(2,1))+  stat_ellipse() +  ggtitle("Weighted Unifrac for Pinacate by years")+themebetaplots+ scale_color_manual(values=colastempbeta)
BetaWpinacate_Chamela
ggsave("news/Weighted Unifrac for betacross by betaplus.png",dpi = 300, units = c("in"), height = 8, width = 10)
BetaUWpinacate_Chamela<-plot_ordination(lepto_sub1_PIN_CHA, ordUW_PIN_CHA, color = "Betaplus",axes=c(2,1))+  stat_ellipse() +  ggtitle("UnWeighted UniFrac for Pinacate by years")+themebetaplots+ scale_color_manual(values=colastempbeta)
BetaUWpinacate_Chamela
ggsave("news/UnWeighted UniFrac for betacross by betaplus.png",dpi = 300, units = c("in"), height = 8, width = 10)

######  UW centroids plot Pinacate-Chamela by years####
metadata_PIN_CHA<-data.frame(sample_data(lepto_sub1_PIN_CHA))
UWyears_PIN_CHA<-data.frame(ordUW_PIN_CHA$vectors[,1],
                            ordUW_PIN_CHA$vectors[,2])
colnames(UWyears_PIN_CHA)[1]<-"MDS1"
colnames(UWyears_PIN_CHA)[2]<-"MDS2"
UWyearscentroid_PIN_CHA<-cbind(UWyears_PIN_CHA,metadata_PIN_CHA)
UWyearscentroid_PIN_CHA
centroids_UWyears_PIN_CHA <- as.data.frame(UWyearscentroid_PIN_CHA %>% 
                                             dplyr::group_by(Betaplus) %>% # calculate functions below for each group
                                             dplyr::summarise(mean_MDS1=mean(MDS1),
                                                              mean_MDS2=mean(MDS2),
                                                              n_MDS1=length(MDS1),
                                                              n_MDS2=length(MDS2),
                                                              stdv_MDS1=sd(MDS1),
                                                              stdv_MDS2=sd(MDS2),
                                                              se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                              se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
UWyears_unifrac_plot_by_betaplus <- ggplot(data = centroids_UWyears_PIN_CHA, aes(x=mean_MDS1, y=mean_MDS2, color=Betaplus))+
  geom_point(size=3) +
  geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                    ymax=mean_MDS2+se_MDS2), width=0.02,  size=0.3)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.02, size=0.3) +
  geom_point(data = UWyearscentroid_PIN_CHA, aes(x=MDS1, y=MDS2, color=Betaplus), alpha=0.2, size=4) +
  scale_color_manual(values = colastempbeta)+
  labs(x = "MDS1 [15.9%]", y="MDS2 [12.8%]")+
  ggtitle("Unweighted UniFrac for Betaplus") +
  themebetaplots
UWyears_unifrac_plot_by_betaplus+stat_ellipse(data=UWyearscentroid_PIN_CHA, aes(x=MDS1, y=MDS2, color=Betaplus),inherit.aes = FALSE)
ggsave("news/UWyears_unifrac_plot_by_betaplus.png",dpi = 300, units = c("in"), height = 8, width = 10)

######  W centroids plot Pinacate-Chamela by years####
Wyears_PIN_CHA<-data.frame(ordW_PIN_CHA$vectors[,1],
                           ordW_PIN_CHA$vectors[,2])
colnames(Wyears_PIN_CHA)[1]<-"MDS1"
colnames(Wyears_PIN_CHA)[2]<-"MDS2"
Wyearscentroid_PIN_CHA<-cbind(Wyears_PIN_CHA,metadata_PIN_CHA)
Wyearscentroid_PIN_CHA
centroids_Wyears_PIN_CHA <- as.data.frame(Wyearscentroid_PIN_CHA %>% 
                                            dplyr::group_by(Betaplus) %>% # calculate functions below for each group
                                            dplyr::summarise(mean_MDS1=mean(MDS1),
                                                             mean_MDS2=mean(MDS2),
                                                             n_MDS1=length(MDS1),
                                                             n_MDS2=length(MDS2),
                                                             stdv_MDS1=sd(MDS1),
                                                             stdv_MDS2=sd(MDS2),
                                                             se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                             se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
Wyears_unifrac_plot_by_betaplus <- ggplot(data = centroids_Wyears_PIN_CHA, aes(x=mean_MDS1, y=mean_MDS2, color=Betaplus))+
  geom_point(size=3) +
  geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                    ymax=mean_MDS2+se_MDS2), width=0.02,  size=0.3)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.02, size=0.3) +
  geom_point(data = Wyearscentroid_PIN_CHA, aes(x=MDS1, y=MDS2, color=Betaplus), alpha=0.2, size=4) +
  scale_color_manual(values = colastempbeta)+
  labs(x = "MDS1 [27.7%]", y="MDS2 [19.3%]")+
  ggtitle("weighted UniFrac for Betaplus") +
  themebetaplots
Wyears_unifrac_plot_by_betaplus+stat_ellipse(data=Wyearscentroid_PIN_CHA, aes(x=MDS1, y=MDS2, color=Betaplus),inherit.aes = FALSE)
ggsave("news/Wyears_unifrac_plot_by_betaplus.png",dpi = 300, units = c("in"), height = 8, width = 10)

###permanovas localities and biomes####
library(vegan)
metadfall<-data.frame(sample_data(lepto_sub1))
metadfall
betacrossUWall<-adonis(DistUW~Locality*year, data=metadfall, permutations = 9999)
betacrossUWall
betacrossWall<-adonis(DistW~Locality*year, data=metadfall, permutations = 9999)
betacrossWall

###perdist2 W localities####

mod00 <- betadisper(DistW, metadfall$Locality, type = "centroid")
mod00
permutest(mod00, permutations = 9999)
anova(mod00)

W_disper_Locality<-boxplot(mod00, col=colas, ylab = "Distance to centroid", xlab="Locality")

###perdist2 UW localities####
mod0 <- betadisper(DistUW, metadfall$Locality, type = "centroid")
mod0
permutest(mod0, permutations = 9999)
anova(mod0)
TukeyHSD(mod0)

UW_disper_Locality<-boxplot(mod0, col=colas, ylab = "Distance to centroid", xlab="Locality")

###permanovas PINACATE CHAMELA YEARS#####

metadata_PIN_CHA<-data.frame(sample_data(lepto_sub1_PIN_CHA))
betayearW_PIN_CHA<-adonis(DistW_PIN_CHA~Betaplus, data=metadata_PIN_CHA,permutations = 9999)
betayearW_PIN_CHA
betayearUW_PIN_CHA<-adonis(DistUW_PIN_CHA~Betaplus, data=metadata_PIN_CHA,permutations = 9999)
betayearUW_PIN_CHA
###perdist2 UW betaplus####
model33 <- betadisper(DistUW_PIN_CHA, metadata_PIN_CHA$Betaplus, type = "centroid")
model33
permutest(model33, permutations = 9999)
anova(model33)
TukeyHSD(model33)

model33

UW_PIN_CHA<-boxplot(model33, col=colastemp, ylab = "UW: Distance to centroid", title(main ="UW Distances to centroid by year"))
UW_PIN_CHA

###perdist2 W betaplus####
model222 <- betadisper(DistW_PIN_CHA, metadata_PIN_CHA$Betaplus, type = "centroid")
model222
permutest(model222, permutations = 9999)
anova(model222)
TukeyHSD(model222)

W_PIN_CHA<-boxplot(model222, col=colastemp, ylab = "W: Distance to centroid", title(main ="W Distances to centroid by year"))

#create ordination
ordW_PIN_CHA = ordinate(lepto_sub1_PIN_CHA, method = "PCoA", distance = DistW_PIN_CHA)
ordUW_PIN_CHA= ordinate(lepto_sub1_PIN_CHA, method = "PCoA", distance = DistUW_PIN_CHA)

plot_scree(ordUW_PIN_CHA)
plot_scree(ordW_PIN_CHA)
ordW_PIN_CHA

BetaWpinacate_Chamela<-plot_ordination(lepto_sub1_PIN_CHA, ordW_PIN_CHA, color = "Betaplus",axes=c(2,1))+  stat_ellipse() +  ggtitle("Weighted Unifrac for Pinacate by years")+themebetaplots+ scale_color_manual(values=colastemp)
BetaWpinacate_Chamela
ggsave("news/Weighted Unifrac for betacross by betaplus.png")
BetaUWpinacate_Chamela<-plot_ordination(lepto_sub1_PIN_CHA, ordUW_PIN_CHA, color = "Betaplus",axes=c(2,1))+  stat_ellipse() +  ggtitle("UnWeighted UniFrac for Pinacate by years")+themebetaplots+ scale_color_manual(values=colastemp)
BetaUWpinacate_Chamela
ggsave("news/UnWeighted UniFrac for betacross by betaplus.png")
