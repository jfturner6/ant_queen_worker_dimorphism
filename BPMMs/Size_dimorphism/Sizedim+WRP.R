### QUEEN-WORKER SIZE DIMORPHISM AND WRP ###

.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))
setwd("/drives/4tb/Juliet/since_Nov_24/regressions/Size_dimorphism/Sizedim+WRP")

library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(geiger)
library(MuMIn)
library(MCMCglmm)
#library(arm)
#library(lme4)
library(purrr)
library(mulTree)
#install.packages("remotes") #Packages needed to install mulTree
#remotes::install_github("TGuillerme/mulTree") #Packages needed to install mulTree
#Further advice for installing mulTree packages can be found on https://github.com/TGuillerme/mulTree

ant_data <- read.csv("merged_data_Jan25.csv", header=T)

# merge queen columns - Queen_HW_mean is from AntWeb, Queen_HW_mean.1 is from elsewhere in literature
queen1 <- dplyr::filter(ant_data,Queen_HW_mean > 0) %>% 
  dplyr::select(animal, Queen_HW_mean) %>% 
  rename(Queen_HW_mean_merged = Queen_HW_mean) 
queen2 <- dplyr::filter(ant_data,Queen_HW_mean.1 > 0)%>%
  dplyr::select(animal, Queen_HW_mean.1) %>%
  rename(Queen_HW_mean_merged = Queen_HW_mean.1)
queen_combined<-rbind(queen1, queen2)
queen_combined_mean<- queen_combined %>%
  group_by(animal) %>%
  mutate(Queen_HW_mean_merged = mean(Queen_HW_mean_merged)) %>%
  distinct(animal, .keep_all = T)
ant_data_with_queen_hw<- left_join(queen_combined_mean, ant_data, by = "animal")

##

#trees from Economo https://datadryad.org/stash/dataset/doi:10.5061/dryad.g579t7k
NCuniform_stem_posterior <- read.tree("15k_NCuniform_stem_posterior.tre")
NCuniform_crown_posterior <- read.tree("15k_NCuniform_crown_posterior.tre")
FBD_stem_posterior <- read.tree("15K_FBD_stem_posterior.tre")
FBD_crown_posterior <- read.tree("15K_FBD_crown_posterior.tre")
all_400_trees <-c(NCuniform_stem_posterior, NCuniform_crown_posterior, FBD_stem_posterior, FBD_crown_posterior)


ant_data<-ant_data_with_queen_hw
ant_data$colony.size <- gsub(",", ".", ant_data$colony.size, fixed = T)
ant_data$colony.size<- as.numeric(ant_data$colony.size)

#### COLONY SIZE ####
## filter (I did not filter out gamergate species here)
## CONFIDENCE >= 0 meaning this analysis will run on all species 
ant_data <- dplyr::filter(ant_data, type == 'ant', Inherent_Constraints_Scale >=0, Min_size >=0, Queen_HW_mean_merged >=0, Supercolonial.Helantera.review_no.stringent.nonstringent <1, social.parasite.Antwiki_no.yes < 1, clonal < 1, hybridisation < 1, Taxonomic_info.Antcat.registered_undescribed.species_Antcat.unregistered_incertae.sedis._subspecies == 0)
ant_data <-dplyr::select(ant_data, animal, Genus, type, Inherent_Constraints_Scale, Min_size, Queen_HW_mean_merged)
ant_data$Genus <- sapply(strsplit(ant_data$Genus, "-"), "[[",1)
ant_data$Max.Reprod.DOL<-c(ant_data$Queen_HW_mean_merged/ant_data$Min_size)
#subset:
ant_data$log10.Max.Reprod.DOL<-log10(ant_data$Max.Reprod.DOL) 

ant_data$animal <- gsub("_", ".", ant_data$animal, fixed = T)
ant_data<-filter(ant_data, animal %in% all_400_trees[[1]]$tip.label)
View(ant_data) # XX species

pruned_all_400_trees<- drop.tip.multiPhylo(all_400_trees, setdiff(all_400_trees[[1]]$tip.label, ant_data$animal)) 


ant_data2 <- ant_data %>% ungroup() %>% as.data.frame() #contained grouping metadata not recognised by multree

############################################################
#MCMCglmm with mulTree
############################################################

#Prepare the data for the mulTree function
mulTree_Reprod_data_sizedim_WRP <- as.mulTree(data = ant_data2, tree = pruned_all_400_trees,
                                             taxa = "animal")

multree_formula_sizedim_WRP <- log10.Max.Reprod.DOL ~ Inherent_Constraints_Scale
# Does worker reproductive potential predict queen-worker size dimorphism?

nitt <- 11000000 
burnin <- 1000000 
thin <- 5000 # parameters based on caste number methods
mul_parameters <- c(nitt, thin, burnin)
set.seed(123)


prior_inv_wishart <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)


mulTree(mulTree.data = mulTree_Reprod_data_sizedim_WRP, formula = multree_formula_sizedim_WRP, priors = prior_inv_wishart,
        parameters = mul_parameters, output = "400_sizedim_WRP", ESS = 1000, convergence = 1.1,
        chains = 2, family = "gaussian")
