### PHYLOGENETIC PATH ANALYSIS ANALYSING WORKER SEX AS A BINARY VARIABLE
# Juliet F. R. Turner, 2025 


## repeated for all 4 MCC trees. 

library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(ggplot2)
library(phylopath)
library(car)
library(dplyr)

ant_data <- read.csv("Nov24_data.csv", header=T) ####### HIGHEST LEVEL CONFIDENCE SUBSET


#fix formatting
ant_data$colony.size<- gsub(",", ".", ant_data$colony.size, fixed = T)
ant_data$colony.size<-as.numeric(ant_data$colony.size)

ant_data$queen.mating.frequency.original<- gsub(",", ".", ant_data$queen.mating.frequency.original, fixed = T)
ant_data$queen.mating.frequency.original<-as.numeric(ant_data$queen.mating.frequency.original)

ant_data$queen.number.continuous<- gsub(",", ".", ant_data$queen.number.continuous, fixed = T)
ant_data$queen.number.continuous<-as.numeric(ant_data$queen.number.continuous)

reprod_path_data <- dplyr::filter(ant_data, type == 'ant', queen.number.continuous >=0, queen.mating.frequency.original >=0, colony.size >=0, 
                                  Inherent_Constraints_Scale >=0, Constraints_Confidence >=0, Supercolonial.Helantera.review_no.stringent.nonstringent <1, social.parasite.Antwiki_no.yes < 1, clonal < 1, hybridisation < 1, Taxonomic_info.Antcat.registered_undescribed.species_Antcat.unregistered_incertae.sedis._subspecies == 0)

reprod_path_data_subset<-dplyr::select(reprod_path_data, animal, type, colony.size, queen.mating.frequency.original, queen.number.continuous, Inherent_Constraints_Scale)
reprod_path_data_subset<-as.data.frame(reprod_path_data_subset)

# MAKE SEX BINARY 
class(reprod_path_data_subset$Inherent_Constraints_Scale) #integer

reprod_path_data_subset$Inherent_Constraints_Scale <- as.numeric(reprod_path_data_subset$Inherent_Constraints_Scale)

reprod_path_data_subset$Inherent_Constraints_Scale <- ifelse(reprod_path_data_subset$Inherent_Constraints_Scale == "1"| 
                                                               reprod_path_data_subset$Inherent_Constraints_Scale == "2", "0", 
                                                             "1") # 0 = SEX
reprod_path_data_subset$totipotency_binary <- as.factor(reprod_path_data_subset$Inherent_Constraints_Scale)

#56 SPECIES

#in worker paper - log10 transformations for CS and MF variables and a square root transformation for worker_polymorphism
reprod_path_data_subset$log10.colony.size<-log10(reprod_path_data_subset$colony.size)
reprod_path_data_subset$log10.eff.mating.freq <-log10(reprod_path_data_subset$queen.mating.frequency.original)
reprod_path_data_subset$log10.queen.number.continuous<-log10(reprod_path_data_subset$queen.number.continuous)
reprod_path_data_subset$animal <- gsub("_", ".", reprod_path_data_subset$animal, fixed = T)
rownames(reprod_path_data_subset) <- reprod_path_data_subset$animal
all_variables<-dplyr::select(reprod_path_data_subset,  animal, type, log10.colony.size, log10.eff.mating.freq, log10.queen.number.continuous, totipotency_binary)
# Convert tibble to a base R data frame
all_variables <- as.data.frame(all_variables)
# Set the row names to the species names
rownames(all_variables) <- all_variables$animal ## 56 species


#Read in phylogenetic trees
NCuniform_stem <- read.tree(file = "15k_NCuniform_stem_mcc.tre")
NCuniform_crown <- read.tree(file = "15K_NCuniform_crown_mcc.tre")
FBD_stem <- read.tree(file = "15K_FBD_stem_mcc.tre")
FBD_crown <- read.tree(file = "15K_FBD_crown_mcc.tre")

#Prune tree
NCuniform_stem_pruned <- drop.tip(NCuniform_stem, setdiff(NCuniform_stem$tip.label, all_variables$animal))
NCuniform_crown_pruned <- drop.tip(NCuniform_crown, setdiff(NCuniform_crown$tip.label, all_variables$animal))
FBD_stem_pruned <- drop.tip(FBD_stem, setdiff(FBD_stem$tip.label, all_variables$animal))
FBD_crown_pruned <- drop.tip(FBD_crown, setdiff(FBD_crown$tip.label, all_variables$animal))

#Rename the variables used in the analysis
all_variables <- all_variables %>% 
  rename(Mating_frequency=log10.eff.mating.freq,
         Colony_size=log10.colony.size,
         Loss_of_Sex=totipotency_binary,
         Queen_Number=log10.queen.number.continuous)

#Prune database
all_variables <- filter(all_variables, animal %in% NCuniform_stem_pruned$tip.label)

#models<- define_model_set(  a = c(Sterility ~ Colony_size, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_Number),  
#                      b = c(Sterility ~ Colony_size, Mating_frequency ~ Colony_size, Sterility ~ Mating_frequency, Mating_frequency ~ Queen_Number),
#                      c =c(Colony_size ~ Sterility, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_Number),
#                      d =c(Sterility ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_Number),
#                      e = c(Mating_frequency ~ Colony_size, Sterility ~ Colony_size, Sterility ~ Queen_Number),
#                      f = c(Mating_frequency ~ Colony_size, Sterility ~ Colony_size, Queen_Number ~ Sterility),
#                      g = c(Mating_frequency ~ Colony_size, Colony_size ~ Sterility, Queen_Number ~ Sterility),
#                      h = c(Mating_frequency ~ Colony_size, Sterility ~ Colony_size, Colony_size ~ Queen_Number),
#                     i = c(Mating_frequency ~ Colony_size, Sterility ~ Colony_size, Queen_Number ~ Colony_size),
#                      j = c(Mating_frequency ~ Colony_size, Colony_size ~ Sterility, Queen_Number ~ Colony_size),
#                      k = c(Mating_frequency ~ Colony_size, Colony_size ~ Queen_Number, Sterility ~ Queen_Number),
#                      l = c(Mating_frequency ~ Colony_size, Colony_size ~ Queen_Number, Queen_Number ~ Sterility),
#                      m = c(Mating_frequency ~ Colony_size, Sterility ~ Colony_size), #REMOVE QN
#                      n = c(Mating_frequency ~ Colony_size, Colony_size ~ Sterility),
#                      o = c(Mating_frequency ~ Colony_size),
#                      p = c(Colony_size ~ Mating_frequency)) 

models<- define_model_set(  a = c(Loss_of_Sex ~ Colony_size, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_Number),  
                            b = c(Colony_size ~ Loss_of_Sex, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_Number),
                            c = c(Loss_of_Sex ~ Colony_size, Loss_of_Sex ~ Queen_Number, Mating_frequency ~ Queen_Number, Mating_frequency ~ Colony_size),
                            d = c(Loss_of_Sex ~ Queen_Number, Mating_frequency ~ Queen_Number, Mating_frequency ~ Colony_size)) 


plot_model_set(models)


##### NCuniform_stem ######
NCuniform_stem_result <- phylo_path(models, data = all_variables, tree = NCuniform_stem_pruned, model = "lambda")
summary(NCuniform_stem_result)
plot(summary(NCuniform_stem_result))
NCuniform_stem_result_average_model_full <- average(NCuniform_stem_result) #, avg_method = "full"

plot(NCuniform_stem_result_average_model_full, algorithm = 'sugiyama', 
     curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, 
     labels = c(Colony_size = "CS", Mating_frequency = "MF", Loss_of_Sex = "S", Queen_Number = "QN"))

NC_stem_plot <- plot(NCuniform_stem_result_average_model_full, algorithm = 'sugiyama', 
                     curvature = 0.04, box_x = 10, box_y = 10, text_size = 5,
                     labels = c(Colony_size = "CS", Mating_frequency = "MF", Loss_of_Sex = "S", Queen_Number = "QN"))

#Anna code for plotting coefficients and CIs
coef_plot(NCuniform_stem_result_average_model_full, reverse_order = TRUE) +
  ggplot2::coord_flip() +
  ggplot2::theme_bw()


##### NCuniform_crown #####
NCuniform_crown_result <- phylo_path(models, data = all_variables, tree = NCuniform_crown_pruned, model = "lambda")
summary(NCuniform_crown_result)
plot(summary(NCuniform_crown_result))
NCuniform_crown_result_average_model_full <- average(NCuniform_crown_result) #, avg_method = "full"

plot(NCuniform_crown_result_average_model_full, algorithm = 'sugiyama', 
     curvature = 0.04, box_x = 10, box_y = 10, text_size = 5,
     labels = c(Colony_size = "CS", Mating_frequency = "MF", Loss_of_Sex = "S", Queen_Number = "QN"))

NC_crown_plot <- plot(NCuniform_crown_result_average_model_full, algorithm = 'sugiyama', 
                      curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, 
                      labels = c(Colony_size = "CS", Mating_frequency = "MF", Loss_of_Sex = "S", Queen_Number = "QN"))




##### FBD_stem #####
FBD_stem_result <- phylo_path(models, data = all_variables, tree = FBD_stem_pruned, model = "lambda")
summary(FBD_stem_result)
plot(summary(FBD_stem_result))
FBD_stem_result_average_model_full <- average(FBD_stem_result, avg_method = "full")

plot(FBD_stem_result_average_model_full, algorithm = 'sugiyama', 
     curvature = 0.04, box_x = 10, box_y = 10, text_size = 5,
     labels = c(Colony_size = "CS", Mating_frequency = "MF", Loss_of_Sex = "S", Queen_Number = "QN"))

FBD_stem_plot <- plot(FBD_stem_result_average_model_full, algorithm = 'sugiyama', 
                      curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, 
                      labels = c(Colony_size = "CS", Mating_frequency = "MF", Loss_of_Sex = "S", Queen_Number = "QN"))

##### FBD_crown #####
FBD_crown_result <- phylo_path(models, data = all_variables, tree = FBD_crown_pruned, model = "lambda")
summary(FBD_crown_result)
plot(summary(FBD_crown_result))
FBD_crown_result_average_model_full <- average(FBD_crown_result) #, avg_method = "full"

plot(FBD_crown_result_average_model_full, algorithm = 'sugiyama', 
     curvature = 0.2, box_x = 10, box_y = 10, text_size = 5, 
     labels = c(Colony_size = "CS", Mating_frequency = "MF", Loss_of_Sex = "S", Queen_Number = "QN"))

FBD_crown_plot <- plot(FBD_crown_result_average_model_full, algorithm = 'sugiyama', 
                       curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, 
                       labels = c(Colony_size = "CS", Mating_frequency = "MF", Loss_of_Sex = "S", Queen_Number = "QN"))

#Anna code for plotting coefficients and CIs
coef_plot(NCuniform_stem_result_average_model_full, reverse_order = TRUE) +
  ggplot2::coord_flip() +
  ggplot2::theme_bw()


# add more options for MF and sterility?
#############################################################



###Summarise the models
##k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight, p=p-value

#NCuniform_stem
NCuniform_stem_summary <- NCuniform_stem_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(NCuniform_stem_summary) <- NULL
NCuniform_stem_summary <- NCuniform_stem_summary %>%
  mutate(phylogeny = "NC uniform stem") %>%
  select(phylogeny, everything())
#write.csv(NCuniform_stem_summary, file = "NCuniform_stem_summary.csv", row.names = FALSE)

#NCuniform_crown
NCuniform_crown_summary <- NCuniform_crown_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(NCuniform_crown_summary) <- NULL
NCuniform_crown_summary <- NCuniform_crown_summary %>%
  mutate(phylogeny = "NC uniform crown") %>%
  select(phylogeny, everything())
# write.csv(NCuniform_crown_summary, file = "NCuniform_crown_summary.csv", row.names = FALSE)

#FBD_stem
FBD_stem_summary <- FBD_stem_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(FBD_stem_summary) <- NULL
FBD_stem_summary <- FBD_stem_summary %>%
  mutate(phylogeny = "FBD stem") %>%
  select(phylogeny, everything())
# write.csv(FBD_stem_summary, file = "FBD_stem_summary.csv", row.names = FALSE)

#FBD_crown
FBD_crown_summary <- FBD_crown_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(FBD_crown_summary) <- NULL
FBD_crown_summary <- FBD_crown_summary %>%
  mutate(phylogeny = "FBD crown") %>%
  select(phylogeny, everything())
# write.csv(FBD_crown_summary, file = "FBD_crown_summary.csv", row.names = FALSE)

##Combine the four data frames
combined_summaries <- rbind(NCuniform_stem_summary, NCuniform_crown_summary, FBD_stem_summary, FBD_crown_summary)
# write.csv(combined_summaries, file = "Path_analysis_summary_siz_var_cont.csv", row.names = FALSE)




######################################################################
  # "_" is "affects", not "affected by"



#Function for extracting Path coefficient, SE, and 95% confidence interval from the Bootstrapped model results
generate_stats <- function(result, value, MCC) {
  CS_S <- c(round(result$coef["Colony_size", "Loss_of_Sex"], 2), 
            round(result$se["Colony_size", "Loss_of_Sex"], 2), 
            round(result$lower["Colony_size", "Loss_of_Sex"], 2), 
            round(result$upper["Colony_size", "Loss_of_Sex"], 2))
  
  S_CS <- c(round(result$coef["Loss_of_Sex", "Colony_size"], 2), 
            round(result$se["Loss_of_Sex", "Colony_size"], 2), 
            round(result$lower["Loss_of_Sex", "Colony_size"], 2), 
            round(result$upper["Loss_of_Sex", "Colony_size"], 2))
  
  CS_MF <- c(round(result$coef["Colony_size", "Mating_frequency"], 2), 
             round(result$se["Colony_size", "Mating_frequency"], 2), 
             round(result$lower["Colony_size", "Mating_frequency"], 2), 
             round(result$upper["Colony_size", "Mating_frequency"], 2))
  
  QN_MF <- c(round(result$coef["Queen_Number", "Mating_frequency"], 2), 
             round(result$se["Queen_Number", "Mating_frequency"], 2), 
             round(result$lower["Queen_Number", "Mating_frequency"], 2), 
             round(result$upper["Queen_Number", "Mating_frequency"], 2))
  
  QN_S <- c(round(result$coef["Queen_Number", "Loss_of_Sex"], 2), 
            round(result$se["Queen_Number", "Loss_of_Sex"], 2), 
            round(result$lower["Queen_Number", "Loss_of_Sex"], 2), 
            round(result$upper["Queen_Number", "Loss_of_Sex"], 2))
  
  stats_noname <- as.data.frame(t(data.frame(CS_affects_S = CS_S, 
                                             S_affects_CS = S_CS,
                                             CS_affects_MF = CS_MF,
                                             QN_affects_MF = QN_MF,
                                             QN_affects_S = QN_S)))
  stats <- rename(stats_noname, 
                  Path_coefficient = V1, 
                  SE = V2, 
                  Lower_95_CI = V3, 
                  Upper_95_CI = V4)
  
  # Add the column "Model_number" to the first position
  stats <- mutate(stats, Model_number := value, Tree := MCC)
  
  # Add row names to a new column
  stats <- rownames_to_column(stats, var = "Causal_relationship")
  
  return(stats)
}






### Calculate confidence intervals ###

### NCUNIFORM STEM ###
NCuniform_stem_result_a <- choice(NCuniform_stem_result, "a", boot = 500) 
NCuniform_stem_result_b <- choice(NCuniform_stem_result, "b", boot = 500)
NCuniform_stem_result_c <- choice(NCuniform_stem_result, "c", boot = 500)
NCuniform_stem_result_d <- choice(NCuniform_stem_result, "d", boot = 500)
NCuniform_stem_result_average_model_full #CI for the average model

NCuniform_stem_coef_stats_a <- generate_stats(result = NCuniform_stem_result_a, value = "a", MCC = "NCuniform_stem")
NCuniform_stem_coef_stats_b <- generate_stats(result = NCuniform_stem_result_b, value = "b", MCC = "NCuniform_stem")
NCuniform_stem_coef_stats_c <- generate_stats(result = NCuniform_stem_result_c, value = "c", MCC = "NCuniform_stem")
NCuniform_stem_coef_stats_d <- generate_stats(result = NCuniform_stem_result_d, value = "d", MCC = "NCuniform_stem")

NCuniform_stem_coef_stats_all_mod <- rbind(NCuniform_stem_coef_stats_a, NCuniform_stem_coef_stats_b, NCuniform_stem_coef_stats_c, 
                                           NCuniform_stem_coef_stats_d)
#######



### NCUNIFORM CROWN ###
NCuniform_crown_result_a <- choice(NCuniform_crown_result, "a", boot = 500) 
NCuniform_crown_result_b <- choice(NCuniform_crown_result, "b", boot = 500)
NCuniform_crown_result_c <- choice(NCuniform_crown_result, "c", boot = 500)
NCuniform_crown_result_d <- choice(NCuniform_crown_result, "d", boot = 500)
NCuniform_crown_result_average_model_full #CI for the average model

NCuniform_crown_coef_stats_a <- generate_stats(result = NCuniform_crown_result_a, value = "a", MCC = "NCuniform_crown")
NCuniform_crown_coef_stats_b <- generate_stats(result = NCuniform_crown_result_b, value = "b", MCC = "NCuniform_crown")
NCuniform_crown_coef_stats_c <- generate_stats(result = NCuniform_crown_result_c, value = "c", MCC = "NCuniform_crown")
NCuniform_crown_coef_stats_d <- generate_stats(result = NCuniform_crown_result_d, value = "d", MCC = "NCuniform_crown")
NCuniform_crown_coef_stats_all_mod <- rbind(NCuniform_crown_coef_stats_a, NCuniform_crown_coef_stats_b, NCuniform_crown_coef_stats_c, 
                                            NCuniform_crown_coef_stats_d)

#######


### FBD CROWN ##
FBD_crown_result_a <- choice(FBD_crown_result, "a", boot = 500) 
FBD_crown_result_b <- choice(FBD_crown_result, "b", boot = 500)
FBD_crown_result_c <- choice(FBD_crown_result, "c", boot = 500)
FBD_crown_result_d <- choice(FBD_crown_result, "d", boot = 500)
FBD_crown_result_average_model_full #CI for the average model

FBD_crown_coef_stats_a <- generate_stats(result = FBD_crown_result_a, value = "a", MCC = "FBD_crown")
FBD_crown_coef_stats_b <- generate_stats(result = FBD_crown_result_b, value = "b", MCC = "FBD_crown")
FBD_crown_coef_stats_c <- generate_stats(result = FBD_crown_result_c, value = "c", MCC = "FBD_crown")
FBD_crown_coef_stats_d <- generate_stats(result = FBD_crown_result_d, value = "d", MCC = "FBD_crown")

FBD_crown_coef_stats_all_mod <- rbind(FBD_crown_coef_stats_a, FBD_crown_coef_stats_b, FBD_crown_coef_stats_c, 
                                      FBD_crown_coef_stats_d)

#######


### FBD STEM ###
FBD_stem_result_a <- choice(FBD_stem_result, "a", boot = 500) 
FBD_stem_result_b <- choice(FBD_stem_result, "b", boot = 500)
FBD_stem_result_c <- choice(FBD_stem_result, "c", boot = 500)
FBD_stem_result_d <- choice(FBD_stem_result, "d", boot = 500)

FBD_stem_result_average_model_full #CI for the average model

FBD_stem_coef_stats_a <- generate_stats(result = FBD_stem_result_a, value = "a", MCC = "FBD_stem")
FBD_stem_coef_stats_b <- generate_stats(result = FBD_stem_result_b, value = "b", MCC = "FBD_stem")
FBD_stem_coef_stats_c <- generate_stats(result = FBD_stem_result_c, value = "c", MCC = "FBD_stem")
FBD_stem_coef_stats_d <- generate_stats(result = FBD_stem_result_d, value = "d", MCC = "FBD_stem")


FBD_stem_coef_stats_all_mod <- rbind(FBD_stem_coef_stats_a, FBD_stem_coef_stats_b, FBD_stem_coef_stats_c, 
                                     FBD_stem_coef_stats_d)



#######



sterility_coef_stats <- rbind(NCuniform_stem_coef_stats_all_mod, NCuniform_crown_coef_stats_all_mod, FBD_crown_coef_stats_all_mod, FBD_stem_coef_stats_all_mod)
write.csv(sterility_coef_stats, file = "250625_allspp_sex_Path_Coefficients_Stats.csv", row.names = F)


















