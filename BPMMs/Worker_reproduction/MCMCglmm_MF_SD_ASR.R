########ASR with MCMCglmm using ASR from corHMM######
##Analysing colony size and worker totipotency (sexual capacity/ no sexual capacity)
#Juliet September 2024, adapted from Louis Bell-Roberts

.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))
setwd("/drives/4tb/Juliet/since_Nov_24/SD_ASRs/Sizedim_MF/MCMCglmm_sizedim_MF")


#Load packages
library(dplyr)
library(ape)
library(phytools)
library(ggplot2)
library(nloptr,
        lib.loc = "/drives/4tb/Juliet/miniconda3/envs/my_r/lib/R/library/")
library(gridExtra,
        lib.loc = "/drives/4tb/Juliet/miniconda3/envs/my_r/lib/R/library/")
library(corHMM)
library(doParallel)
library(MCMCglmm)

#Run in parallel 
registerDoParallel(40) # 


################
#Custom functions     adjust for 2rat_ARD sterility (1 = no sexual capacity/sterile, 0 = sex/non-sterile)
################

#Calculate the most probable state for each node in the corHMM ASR
##Define function to find index of maximum value in a vector using the max_column function
max_column <- function(x) {
  return(names(x)[which.max(x)])
}

# Define a function to assign sterile/non-sterile to each node in my corHMM ASR - this function is for when using corHMM ASRs with 2 rate categories
assign_node_state <- function(df) {
  df$node_state <- ifelse(df$MaxColumn %in% c("X.1.R1."), "low", "high") 
  return(df)
}


 
# Define function to rename columns as poly2 and animal
rename_cols <- function(df) {
  names(df)[1] <- "toti2"
  names(df)[2] <- "animal"
  return(df)
}


################################################
#Data preparation

#Read in data file
new.data <- read.csv("MF_sizedim_ASR_df.csv") # subset includes both WRP and colony size

ant_data <- new.data

#Read in sample of 400 phylogenetic trees
ant_trees_pruned <- read.tree(file ="MF_sizedim_ASR_tree.tre")

#Add node labels as they were not originally included with the phylo object
##Edges and node labels used later to work out where transitions in sterility happens
for (i in 1:length(ant_trees_pruned)) {
  ant_trees_pruned[[i]]$node.label <- paste("Node", 1:ant_trees_pruned[[i]]$Nnode, sep = "")
}


##########################################################################################
#ANCESTRAL STATE RECONSTRUCTION WITH COLONY SIZE AND TOTIPOTENCY
##########################################################################################

###############
#Read in the 400 corHMM model results
##corHMM files must be read into R in a particular order - the particular tree used for each corHMM model must correspond to the tree used for each MCMCglmm model later on
###Create an empty list to hold the results
corHMM_results <- list()

#Set folder path where corHMM model results are saved - reconstruction of WRP nodes
folder_path <- file.path("/drives/4tb/Juliet/since_Nov_24/SD_ASRs/Sizedim_MF/corHMM_sizedim_MF/sizedim_MF_output")

#Get the file names
file_names <- list.files(folder_path)

#Extract the numbers from the file names using regular expressions
model_numbers <- as.numeric(gsub(".*_(\\d+)\\.rda", "\\1", file_names))

#Sort the file names based on the corresponding numbers
sorted_file_names <- file_names[order(model_numbers)]

#Read in the files and save them to a list
corHMM_results <- lapply(sorted_file_names, function(x) readRDS(file.path(folder_path, x)))

#Extracting the states component from each model and create a new list of just the states
corHMM_states <- lapply(corHMM_results, function(model) model$states)


###############
#For each tree, estimate the most likely state of totipotency (sex/ no sex) for each node of the phylogeny

#First, convert each of the "matrix" "array" objects that are in the corHMM_states list to data frame objects before running the subsequent functions
##Apply function to each matrix in the list
corHMM_states_conv <- lapply(corHMM_states, data.frame)

#Apply max_column custom function to each row of the data frame to find the most likely state at each node in the corHMM ASRs
##Apply function to each data frame in the list
corHMM_states_conv <- lapply(corHMM_states_conv, function(df) {
  df$MaxColumn <- apply(df, 1, max_column)
  return(df)
})

###############

##Apply the function to each dataframe in the list
corHMM_states_conv <- lapply(corHMM_states_conv, assign_node_state)


###############
#Data preparation
##Create "animal.Node" number column pasted before the node number for each row
for (i in seq_along(corHMM_states_conv)) {
  corHMM_states_conv[[i]]$node_number <- paste0("animal.Node", seq_len(nrow(corHMM_states_conv[[i]])))
}

#Replace "animal.Node1" with "(Intercept)" in each data frame in the list
for (i in seq_along(corHMM_states_conv)) {
  corHMM_states_conv[[i]]$node_number <- gsub("\\banimal\\.Node1\\b", "(Intercept)", corHMM_states_conv[[i]]$node_number)
}

#Create new list of data frames extracting just the node_state and node_number columns of 'corHMM_states_conv'
corHMM_states_conv_minimal <- lapply(corHMM_states_conv, function(df) df[, c("node_state", "node_number")])

#Set row names to match the values in node_number column
corHMM_states_conv_minimal <- lapply(corHMM_states_conv_minimal, function(df) {
  rownames(df) <- df[, "node_number"]
  return(df)
})

#Get data for each species tip and its totipotency data
##Create 'nodecode' column" for sxData by pasting animal and then also the species tip
ant_data$nodecode <- paste("animal", ant_data$animal,sep=".")

#Create new data frame extracting just the totipotency_binary (prev. CasteBin) and nodecode columns
ant_data_minimal <- ant_data[, c("dimorphism.binary", "nodecode")]

#Make the row names equal to nodecode column
rownames(ant_data_minimal) <- ant_data_minimal[, "nodecode"]

#The column names of corHMM_states_conv_minimal and sxData_minimal must match for rbind to work
##Use lapply to change column names for each data frame in the list of corHMM_states_conv_minimal using the rename_cols custom function
corHMM_states_conv_minimal <- lapply(corHMM_states_conv_minimal, rename_cols)

#Change column names using names function for sxData_minimal
names(ant_data_minimal) <- c("toti2", "animal")

# ADDED
ant_data_minimal<- ant_data_minimal %>%
  mutate(toti2 = ifelse(toti2 == 0, "low", "high"))

#Append the rows from the 'corHMM_states_conv_minimal' data frame to the 'sxData_minimal' data frame using the rbind operation
##Use lapply to bind sxData_minimal to each element in corHMM_states_conv_minimal
Totipotency <- lapply(corHMM_states_conv_minimal, function(df) {
  rbind(df, ant_data_minimal)
})


#############
#Create Transition Dataset
#Extract edge information for each tree in the multiPhylo object and create a list of dataframes containing information on the relationship between parent and offspring nodes across the tree
TransDat <- lapply(ant_trees_pruned, function(x) as.data.frame(x$edge)) ##Each value appears twice in column V1 as nodes have two descendants

#Subtract tips and create column 'V1name' for each data frame in the TransDat list to identify node names for each row (V1 = parent nodes)
for (i in seq_along(TransDat)) {
  TransDat[[i]]$V1name <- paste("animal.Node", TransDat[[i]]$V1-length(ant_trees_pruned[[1]]$tip.label), sep = "")
}

#Repeat for V2 (offspring nodes/tips)
for (i in seq_along(TransDat)) {
  TransDat[[i]]$V2name <- paste("animal.Node", TransDat[[i]]$V2-length(ant_trees_pruned[[1]]$tip.label), sep = "")
} #'V2name' will include some negative numbers since some of these are tips rather than nodes. Swap these for the correct tip labels later on

#Assign tip number to each species and paste 'animal' in front of it to match with row names in the 'Poly' dataframe created above
treesp <- data.frame(tip.label=paste("animal.", ant_trees_pruned[[1]]$tip.label,sep=""), no=1:length(ant_trees_pruned[[1]]$tip.label)) #While tree topology changes across the sample of 400 trees, the tip number that a species is labelled with always stays the same even if its position on the tree moves

#Replace negative values with tip labels by matching V2 in TransDat with 'no' in treesp - creates the 'descendants' column which contains this information
TransDat_neg_rep <- list()
for (i in seq_along(TransDat)) {
  TransDat_neg_rep[[i]] <- data.frame(TransDat[[i]], ancestors=NA, descendants=treesp$tip.label[match(TransDat[[i]]$V2, treesp$no)])
}

#Assign class 'character' to 'ancestors' and 'descendants'
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {df$ancestors <- as.character(df$ancestors); return(df)})
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {df$descendants <- as.character(df$descendants); return(df)})

#Values for ancestors are currently NAs, so assign V1name 
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {
  df$ancestors <- as.character(ifelse(!is.na(df$ancestors), df$ancestors, df$V1name))
  return(df)
})

#Replace any descendants (i.e. those that aren't tips) with V2 name
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {
  df$descendants <- as.character(ifelse(!is.na(df$descendants), df$descendants, df$V2name))
  return(df)
}) #The names of each ancestor and descendant are now identified. Since they have 'animal.' pasted in front, 
#they will match with row names in the 'Totipotency' (previously 'Poly') dataframe. 


#############
#Combine colony size estimates with transition dataset
##Match up ancestor and descendant colony size predictions based on row names
TransDat_neg_rep <- lapply(seq_along(TransDat_neg_rep), function(i) {
  data.frame(TransDat_neg_rep[[i]], 
             ancTotipotency=Totipotency[[i]]$toti2[match(TransDat_neg_rep[[i]]$ancestors,Totipotency[[i]]$animal)],
             desTotipotency=Totipotency[[i]]$toti2[match(TransDat_neg_rep[[i]]$descendants,Totipotency[[i]]$animal)])
})

#############
#Calculate the number of different types of transitions between caste number
##Make a table of the number of types of descendant each ancestor has
###Possible transition types:
# 1 polymorphic and 1 monomorphic descendant, 2 poly vs. 0 mono, 0 poly vs. 2 mono
obs_list <- lapply(TransDat_neg_rep, function(df) data.frame(table(df$desTotipotency, df$ancestors)))
obs_list <- lapply(obs_list, function(df){
  df$CAT <- df$Freq 
  return(df)
})

#Assign transition category types (CAT) based on the number of castes that each descendant has
for(i in seq_along(obs_list)) {
  obs <- obs_list[[i]]
  
#Var1 is the column of descendant state, Var2 = ancestral node label
obs$CAT[obs$Freq == 2 & obs$Var1 == "high"] <- "only.high"      
obs$CAT[obs$Freq == 0 & obs$Var1 == "low"] <- "only.high"
# nodes with 2 non.sterile descendants:
obs$CAT[obs$Freq == 2 & obs$Var1 == "low"] <- "only.low"
obs$CAT[obs$Freq == 0 & obs$Var1 == "high"] <- "only.low"
# nodes with 1 sterile and 1 non.sterile descendant:
obs$CAT[obs$Freq == 1] <- "both"
  
  # Assign the modified data frame back to the list
  obs_list[[i]] <- obs
}


#############
#Assign Types of Transitions Between sex/no sex to each Node

#Remove NAs in the transitions dataset
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(x) {
  x <- x[!is.na(x$ancTotipotency),] #This deletes node 1 which doesn't have an ancestor
  x <- x[!is.na(x$desTotipotency),] #Deletes rows where totipotency is uncertain
  return(x)
})

#Assign transition types in TransDat based on the values calculated in obs_list
for (i in seq_along(TransDat_neg_rep)) {
  # match ancestors with corresponding CAT value from obs_list
  TransDat_neg_rep[[i]] <- data.frame(TransDat_neg_rep[[i]], CAT = obs_list[[i]]$CAT[match(TransDat_neg_rep[[i]]$ancestors, obs_list[[i]]$Var2)])
}

#Create a new column that classifies nodes based on ancestor state & transition type. 
for (i in seq_along(TransDat_neg_rep)) {
  # create a new column in each data frame and concatenate ancTotipotency and CAT columns
  TransDat_neg_rep[[i]]$CAT2 <- paste(TransDat_neg_rep[[i]]$ancTotipotency, TransDat_neg_rep[[i]]$CAT, sep = ".")
}


#### NEW VERSION
#Assign 'double transitions' (e.g. monomorphic ancestor to only polymorphic descendants) as only single transition events 
for (i in seq_along(TransDat_neg_rep)) {
  # Check if the current element is a data frame
  if (is.data.frame(TransDat_neg_rep[[i]])) {
    # Proceed to recode if it's a data frame
    TransDat_neg_rep[[i]]$CAT2[TransDat_neg_rep[[i]]$CAT2 == "low.only.high"] <- "low.both" 
    TransDat_neg_rep[[i]]$CAT2[TransDat_neg_rep[[i]]$CAT2 == "high.only.low"] <- "high.both"             
  } else {
    # Print a warning message for non-data frame elements
    warning(paste("Element", i, "is not a data frame."))
  }
}

# After the loop, check unique values from the first data frame to see if the changes took effect
print(unique(TransDat_neg_rep[[1]]$CAT2)) ### should now be only four
table(TransDat_neg_rep[[1]]$CAT2) 

                   
                   

###### CHANGED 29/09/2024 this removes last element which may be an empty character ##### CHECK 
#TransDat_neg_rep <- TransDat_neg_rep[-length(TransDat_neg_rep)]

#############
#Add estimates for the continuous trait for the extant species to the transition datasets

#Match up tips with colony size data based on nodecode column
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(x) {
  #print(class(x))
  x <- data.frame(x, continuous.trait = ant_data$queen.mating.frequency.original[match(x$descendants, ant_data$nodecode)])   
  return(x)
})       


#############
#Assign ancestral nodes to the 'animal' column
#When running the model, this will allow us to estimate colony size in 'ancestors'
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {
  df$animal <- gsub("animal.", "", df$ancestors)
  return(df)
})

#Assign transition types as a factor
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {
  df$CAT2 <- as.factor(df$CAT2)
  return(df)
})

#Remove rownames from each data frame in the list
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(x) {
  rownames(x) <- NULL
  return(x)
})

#############
#Model to estimate colony size values across each transition type

#Prior for model with Gaussian response variable
prior1 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))

# Parallel models with outputs saved as .rds files. Remove the global intercept and run over 400 trees. #gaussian assumes normality of response (colony size)
foreach(i = 1:400) %dopar% {
  #1st model
  model1 <- MCMCglmm(continuous.trait ~ CAT2-1, random = ~animal, family = "gaussian", nodes = "ALL", prior = prior1, pedigree = ant_trees_pruned[[i]], data = TransDat_neg_rep[[i]], nitt = 1100000, burnin = 100000, thin = 1000, pr=TRUE, verbose = F)
  #Run 2nd model
  model2 <- MCMCglmm(continuous.trait ~ CAT2-1, random = ~animal, family = "gaussian", nodes = "ALL", prior = prior1, pedigree = ant_trees_pruned[[i]], data = TransDat_neg_rep[[i]], nitt = 1100000, burnin = 100000, thin = 1000, pr=TRUE, verbose = F)
  
  
  # Save the model as an .rds file for 1st and 2nd chain
  saveRDS(model1, file = file.path("1st_run", paste0("MF_SD_1M_100k_1k_", i, "_1stRun.rds")))
  saveRDS(model2, file = file.path("2nd_run", paste0("MF_SD_1M_100k_1k_", i, "_2ndRun.rds")))
  
}

