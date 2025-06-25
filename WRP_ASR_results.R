########################Plot transitions of the corHMM ASR for WRP and calculate summary statistics for each transition type######################
#Juliet October 2024, adapted from Louis Bell-Roberts


# Use conda environment library path
.libPaths("/drives/4tb/Juliet/miniconda3/envs/my_r/lib/R/library") # added 19/03/2025

#Load packages
library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(geiger)
library(corHMM)
library(phangorn)
library(gtools)
library(ggtree)

setwd("/drives/4tb/Juliet/since_Nov_24/WRP_ASRs/WRP_4level_ASR")

#Load custom functions

#Function that takes a matrix or data frame as input and returns a vector of the column indices corresponding to the maximum values in each row - replicates the behaviour of the max.col function
my_max_col <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("Input must be a matrix or data frame")
  }
  apply(x, 1, which.max)
}

#The extract_frequency function is designed to extract frequencies from a list of data frames (table_list) based on a specified pattern found in the 'Var1' column of each data frame. Here's a summary of what the function does:
extract_frequency <- function(table_list, pattern) {
  lapply(table_list, function(df) {
    if (pattern %in% df$Var1) {
      subset_df <- df[grep(pattern, df$Var1), ]
      # Extract the 'Freq' column from the subsetted data frame
      freq <- subset_df$Freq
    } else {
      # If 'pattern' is not present, return 0
      freq <- 0
    }
    return(freq)
  })
}

###############################################################

#Read in data file
sxData <- read.csv("/drives/4tb/Juliet/since_Nov_24/WRP_ASRs/WRP_4level_ASR/ASR_WRP_df.csv")
# 153 species with highest confidence data, only WRP info and species names
sxData <- sxData %>%
  rename(WRP = Inherent_Constraints_Scale)


ant_trees_pruned <- read.tree(file ="/drives/4tb/Juliet/since_Nov_24/WRP_ASRs/WRP_4level_ASR/ASR_WRP_pruned_tree.tre")

###############
###Read in the 400 corHMM model results in ascending order
###############

# create an empty list to hold the results
ASR_models <- list()

# set folder path
folder_path <- ("/drives/4tb/Juliet/since_Nov_24/WRP_ASRs/WRP_4level_ASR/rat2_ARD")

# get the file names
file_names <- list.files(folder_path)

# sort the file names based on the corresponding numbers
sorted_file_names <- mixedsort(file_names)

# read in the files and save them to a list
ASR_models <- lapply(sorted_file_names, function(x) readRDS(file.path(folder_path, x)))

####################################################################################################
############ NEW CODE 24TH FEB 25 -- ROOT STATE PROBABILTIES
# Check first few models
#lapply(ASR_models[1:5], function(model) str(model$root.p)) #would work if not using 'yang' method

############# collapsing into 1 rate category

# Extract and collapse root state probabilities based on two rate categories (rate_cat = 2)
root_probs_list <- lapply(ASR_models, function(model) {
  if (!is.null(model$states)) {
    # Extract root state probabilities for both rate categories (first row)
    rate_1_probs <- model$states[1, 1:4]  # Probabilities for the first rate category
    rate_2_probs <- model$states[1, 5:8]  # Probabilities for the second rate category
    
    # Collapse them into 4 categories by choosing the more likely probability for each state
    collapsed_probs <- pmax(rate_1_probs, rate_2_probs)
    
    # Return the collapsed probabilities (numeric vector of 4 categories)
    return(collapsed_probs)
  } else {
    return(rep(NA, 4))  # Handle missing cases
  }
})

# probability of each root state for each of 400 models
str(root_probs_list)

# Convert the list of root probabilities to a matrix where each row corresponds to a model,
# and each column corresponds to one of the four states (1, 2, 3, 4).
root_probs_matrix <- do.call(rbind, root_probs_list)

# Calculate the average probability for each state (1, 2, 3, 4) across all models.
overall_root_probs <- colMeans(root_probs_matrix)

# Print the overall root probabilities
overall_root_probs

# (1,R1)         (2,R1)       (3,R1)     (4,R1) 
# 0.5904840.   0.1160515.   0.1082929.   0.1280553 

# Combine all root probabilities for all 400 models
combined_probs <- do.call(rbind, root_probs_list)

# Calculate the 95% credible interval for each state (1 to 4)
credible_intervals <- apply(combined_probs, 2, function(x) {
  # Calculate the 2.5th and 97.5th percentiles for the root state probabilities
  quantile(x, probs = c(0.025, 0.975))
})

# Print the credible intervals for each root state
print(credible_intervals)


# PLOTTING

# Load necessary library
library(ggplot2)

# Prepare data for plotting with modified state labels
root_states <- c("Full", "Reduced", "Males only", "No offspring")  # Updated labels
estimates <- overall_root_probs
lower_ci <- c(0.0007061208, 4.194245e-05, 5.954294e-05, 0.002444829)
upper_ci <- c(0.9723850063, 0.6260647114, 0.4535402898, 0.708239963)

# Ensure root_states are treated as a factor with specific order
root_states <- factor(root_states, levels = c("Full", "Reduced", "Males only", "No offspring"))

# Create a data frame for plotting
plot_data <- data.frame(
  RootState = root_states,
  Estimate = estimates,
  LowerCI = lower_ci,
  UpperCI = upper_ci
)

# Plot using ggplot2
ggplot(plot_data, aes(x = RootState, y = Estimate)) +
  geom_point(size = 4, color = "#D33682") +  # Plot points for the estimates
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +  # Add error bars for CI
  labs(
    x = "Root State",
    y = "Estimate"
  ) +
  theme_minimal() +  # Apply a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for clarity
    panel.grid = element_blank(),  # Remove grid lines
    axis.title = element_text(size = 14)  # Increase axis title font size
  )






####################################################################################################

###############
###For each tree, estimate the most likely number of X for each node of the phylogeny
###############

#Convert ASR estimates to list of data frames
ASR_states_list <- lapply(ASR_models, function(df) as.data.frame(df$states))

# Create an empty list to store the results
TransDat_neg_rep_list <- vector(mode = "list", length = length(ASR_states_list))

##############
for (i in seq_along(ASR_states_list)) {
  ASR_states <- ASR_states_list[[i]]
  
  #Combine the probability scores for the two different rate classes for each caste number
  ASR_states_summed <- data.frame(WRP_1 = ASR_states$`(1,R1)` + ASR_states$`(1,R2)`, WRP_2 = ASR_states$`(2,R1)` + ASR_states$`(2,R2)`, WRP_3 = ASR_states$`(3,R1)` + ASR_states$`(3,R2)`, WRP_4 = ASR_states$`(4,R1)` + ASR_states$`(4,R2)`)
  
  #Identify which WRP is the highest probability for each node
  ASR_states_summed$max_col <- colnames(ASR_states_summed)[my_max_col(ASR_states_summed)] #originally was using the max.col() function however it was making rounding errors
  
  #Assign node number to each row of the data frame, given that node number corresponds to row number in the ASR_states_summed dataframe
  ASR_states_summed$node_number <- paste0("animal.Node", seq_len(nrow(ASR_states_summed)))
  
  #Create new list of data frames with just the node_state and node_number columns
  ASR_states_summed <- ASR_states_summed %>% dplyr::select(max_col, node_number)
  
  #######
  
  #Get data for each species tip and its WRP data
  ##Create 'nodecode' column" which pastes animal and then also the species tip
  sxData$nodecode <- paste("animal", sxData$animal, sep=".")
  # Create new data frame with just the WRP and nodecode columns
  sxData_minimal <- sxData %>% dplyr::select(WRP, nodecode)
  #Make the row names equal to nodecode column
  rownames(sxData_minimal) <- sxData_minimal[, "nodecode"]
  
  #Join 'ASR_states_summed' on top of 'sxData_minimal' using rbind
  ##Use lapply to bind sxData_minimal to each element in ASR_states_summed
  ASR_states_summed <- ASR_states_summed %>% dplyr::rename(poly2 = max_col, animal = node_number)
  sxData_minimal <- sxData_minimal %>% dplyr::rename(poly2 = WRP, animal = nodecode)
  
  sxData_minimal$poly2 <- ifelse(sxData_minimal$poly2 == 1, "WRP_1",
                                 ifelse(sxData_minimal$poly2 == 2, "WRP_2",
                                        ifelse(sxData_minimal$poly2 == 3, "WRP_3", #based on dataset, not model output
                                               ifelse(sxData_minimal$poly2 == 4, "WRP_4", NA))))
  
  Poly <- rbind(ASR_states_summed, sxData_minimal)
  
  
  #Turn tree edges into a data frame - displays the parent node (V1) and the two offspring nodes (V2). 653 is the root node number as there are 652 tips in the tree
  TransDat <- as.data.frame(ASR_models[[i]]$phy$edge)
  
  #Name each of the nodes in the phylogeny with correct name
  TransDat$V1name <- paste("animal.Node", TransDat$V1-length(ASR_models[[i]]$phy$tip.label), sep = "")
  
  #Do the same but for the descendant nodes
  TransDat$V2name <- paste("animal.Node", TransDat$V2-length(ASR_models[[i]]$phy$tip.label), sep = "")
  
  #Number each of the species in the tree by their tip number
  treesp <- data.frame(tip.label=paste("animal.", ASR_models[[i]]$phy$tip.label,sep=""), no=1:length(ASR_models[[i]]$phy$tip.label))
  
  #Match species names with the descendant nodes
  TransDat_neg_rep <- data.frame(TransDat, ancestors=NA, descendents=treesp$tip.label[match(TransDat$V2, treesp$no)])
  
  # make sure these are of class 'character'
  TransDat_neg_rep$ancestors <- as.character(TransDat_neg_rep$ancestors)
  TransDat_neg_rep$descendents <- as.character(TransDat_neg_rep$descendents)
  
  # ancestors currently NAs, so give these V1name 
  TransDat_neg_rep$ancestors <- as.character(ifelse(!is.na(TransDat_neg_rep$ancestors),TransDat_neg_rep$ancestors,TransDat_neg_rep$V1name))
  
  # replace any descendents (i.e. those that aren't tips) with V2 name
  TransDat_neg_rep$descendents <- as.character(ifelse(!is.na(TransDat_neg_rep$descendents),TransDat_neg_rep$descendents,TransDat_neg_rep$V2name))
  
  
  ################
  #Step that assigns WRP to ancestor and descendant
  TransDat_neg_rep <- data.frame(TransDat_neg_rep, 
                                 ancPoly=Poly$poly2[match(TransDat_neg_rep$ancestors, Poly$animal)],
                                 desPoly=Poly$poly2[match(TransDat_neg_rep$descendents, Poly$animal)])
  
  #create new column with concatenated values from ancPoly and desPoly
  TransDat_neg_rep$trans <- paste0(TransDat_neg_rep$ancPoly, "_to_", TransDat_neg_rep$desPoly)
  TransDat_neg_rep$trans <- as.factor(TransDat_neg_rep$trans)
  
  #Append TransDat_neg_rep to a list
  TransDat_neg_rep_list[[i]] <- TransDat_neg_rep
}
#######


#############
#Calculate summary statistics for the number of transitions of each type over the 400 trees
#############

#For a single dataframe
table(TransDat_neg_rep_list[[1]]$trans)

#For all 400 data frames

# First, apply the table() function to the 'trans' column of each data frame
table_list <- lapply(TransDat_neg_rep_list, function(df) data.frame(table(df$trans)))

##WRP_1_to_WRP_1
#Calculate the mean frequency of occurrences of 'WRP_1_to_WRP_1' across the table_list list of dataframes. 
# If the column 'WRP_1_to_WRP_1' exists in a data frame, extract the frequencies associated with it and compute the mean. 
#If the column doesn't exist in a data frame, assign a frequency of 0 for that data frame.

freq_WRP_1_to_WRP_1 <- extract_frequency(table_list, "WRP_1_to_WRP_1")
mean(unlist(freq_WRP_1_to_WRP_1)) #sorted tree = 24.3725     # unsorted tree = 61.8975    #highconf = 50.1575       #new Nov = 67
median(unlist(freq_WRP_1_to_WRP_1)) 
range(unlist(freq_WRP_1_to_WRP_1)) 

##WRP_1_to_WRP_2
freq_WRP_1_to_WRP_2 <- extract_frequency(table_list, "WRP_1_to_WRP_2")
mean(unlist(freq_WRP_1_to_WRP_2)) 
median(unlist(freq_WRP_1_to_WRP_2)) 
range(unlist(freq_WRP_1_to_WRP_2)) 

##WRP_1_to_WRP_3
freq_WRP_1_to_WRP_3 <- extract_frequency(table_list, "WRP_1_to_WRP_3")
mean(unlist(freq_WRP_1_to_WRP_3)) 
median(unlist(freq_WRP_1_to_WRP_3)) 
range(unlist(freq_WRP_1_to_WRP_3))

##WRP_1_to_WRP_4
freq_WRP_1_to_WRP_4 <- extract_frequency(table_list, "WRP_1_to_WRP_4")
mean(unlist(freq_WRP_1_to_WRP_4)) 
median(unlist(freq_WRP_1_to_WRP_4)) 
range(unlist(freq_WRP_1_to_WRP_4)) 

##WRP_2_to_WRP_2
freq_WRP_2_to_WRP_2 <- extract_frequency(table_list, "WRP_2_to_WRP_2")
mean(unlist(freq_WRP_2_to_WRP_2)) 
median(unlist(freq_WRP_2_to_WRP_2)) 
range(unlist(freq_WRP_2_to_WRP_2)) 

##WRP_2_to_WRP_3
freq_WRP_2_to_WRP_3 <- extract_frequency(table_list, "WRP_2_to_WRP_3")
mean(unlist(freq_WRP_2_to_WRP_3)) 
median(unlist(freq_WRP_2_to_WRP_3)) 
range(unlist(freq_WRP_2_to_WRP_3)) 

##WRP_2_to_WRP_4
freq_WRP_2_to_WRP_4 <- extract_frequency(table_list, "WRP_2_to_WRP_4") 
mean(unlist(freq_WRP_2_to_WRP_4)) 
median(unlist(freq_WRP_2_to_WRP_4)) 
range(unlist(freq_WRP_2_to_WRP_4)) 

##WRP_3_to_WRP_3
freq_WRP_3_to_WRP_3 <- extract_frequency(table_list, "WRP_3_to_WRP_3")
mean(unlist(freq_WRP_3_to_WRP_3)) 
median(unlist(freq_WRP_3_to_WRP_3)) 
range(unlist(freq_WRP_3_to_WRP_3)) 

##WRP_3_to_WRP_4
freq_WRP_3_to_WRP_4 <- extract_frequency(table_list, "WRP_3_to_WRP_4")
mean(unlist(freq_WRP_3_to_WRP_4)) 
median(unlist(freq_WRP_3_to_WRP_4)) 
range(unlist(freq_WRP_3_to_WRP_4)) 

##WRP_4_to_WRP_4
freq_WRP_4_to_WRP_4 <- extract_frequency(table_list, "WRP_4_to_WRP_4")
mean(unlist(freq_WRP_4_to_WRP_4)) 
median(unlist(freq_WRP_4_to_WRP_4)) 
range(unlist(freq_WRP_4_to_WRP_4)) 

##Losses of WRP

##WRP_2_to_WRP_1
freq_WRP_2_to_WRP_1 <- extract_frequency(table_list, "WRP_2_to_WRP_1")
mean(unlist(freq_WRP_2_to_WRP_1)) 
median(unlist(freq_WRP_2_to_WRP_1)) 
range(unlist(freq_WRP_2_to_WRP_1)) 

##WRP_3_to_WRP_1
freq_WRP_3_to_WRP_1 <- extract_frequency(table_list, "WRP_3_to_WRP_1")
mean(unlist(freq_WRP_3_to_WRP_1)) 
median(unlist(freq_WRP_3_to_WRP_1)) 
range(unlist(freq_WRP_3_to_WRP_1)) 

##WRP_3_to_WRP_2
freq_WRP_3_to_WRP_2 <- extract_frequency(table_list, "WRP_3_to_WRP_2")
mean(unlist(freq_WRP_3_to_WRP_2))
median(unlist(freq_WRP_3_to_WRP_2)) 
range(unlist(freq_WRP_3_to_WRP_2)) 

##WRP_4_to_WRP_1
freq_WRP_4_to_WRP_1 <- extract_frequency(table_list, "WRP_4_to_WRP_1")
mean(unlist(freq_WRP_4_to_WRP_1)) 
median(unlist(freq_WRP_4_to_WRP_1)) 
range(unlist(freq_WRP_4_to_WRP_1)) 

##WRP_4_to_WRP_2
freq_WRP_4_to_WRP_2 <- extract_frequency(table_list, "WRP_4_to_WRP_2")
mean(unlist(freq_WRP_4_to_WRP_2)) 
median(unlist(freq_WRP_4_to_WRP_2)) 
range(unlist(freq_WRP_4_to_WRP_2))

##WRP_4_to_WRP_3
freq_WRP_4_to_WRP_3 <- extract_frequency(table_list, "WRP_4_to_WRP_3")
mean(unlist(freq_WRP_4_to_WRP_3)) 
median(unlist(freq_WRP_4_to_WRP_3)) 
range(unlist(freq_WRP_4_to_WRP_3)) 

##############
## NEW CODE ADDED 16/01/2024
# Combine transition frequencies into a summary dataframe
# Extract unique transition types across all tables
all_transition_types <- unique(unlist(lapply(table_list, function(df) as.character(df$Var1))))

# Create an empty dataframe to store summary statistics
summary_stats <- data.frame(
  Transition = all_transition_types,
  Mean = numeric(length(all_transition_types)),
  Median = numeric(length(all_transition_types)),
  Min = numeric(length(all_transition_types)),
  Max = numeric(length(all_transition_types)),
  stringsAsFactors = FALSE
)

# Calculate statistics for each transition type
for (i in seq_along(all_transition_types)) {
  transition <- all_transition_types[i]
  # Extract frequencies for the current transition type
  freqs <- extract_frequency(table_list, transition)
  freqs_unlist <- unlist(freqs)
  
  # Populate the summary dataframe
  summary_stats[i, "Mean"] <- mean(freqs_unlist)
  summary_stats[i, "Median"] <- median(freqs_unlist)
  summary_stats[i, "Min"] <- min(freqs_unlist)
  summary_stats[i, "Max"] <- max(freqs_unlist)
}

# View the summary dataframe
View(summary_stats)

# Optionally, save the summary dataframe to a CSV file
#write.csv(summary_stats, "transition_summary_statistics.csv", row.names = FALSE)


############# MODES

# Load necessary library
library(dplyr)

# Function to calculate the mode
calculate_mode <- function(x) {
  uniq_vals <- unique(x)
  uniq_vals[which.max(tabulate(match(x, uniq_vals)))]
}

# Extract unique transition types from the first data frame
unique_transitions <- unique(TransDat_neg_rep_list[[i]]$trans) #### modified 17JAN25

# Initialize an empty list to store the mode for each transition type
mode_frequencies <- list()

# Loop through each unique transition type, calculate the mode, and store it
for (transition in unique_transitions) {
  # Extract frequencies for the current transition type across all data frames
  freq_current_transition <- extract_frequency(table_list, transition)
  
  # Flatten the list into a single vector
  freq_vector <- unlist(freq_current_transition)
  
  # Calculate the mode for the current transition type
  mode_freq <- calculate_mode(freq_vector)
  
  # Store the result in the list
  mode_frequencies[[transition]] <- mode_freq
}

# Convert the list of mode frequencies to a dataframe
transition_summary_df <- data.frame(
  Transition = names(mode_frequencies),
  ModeFrequency = unlist(mode_frequencies)
)

# View the resulting dataframe
View(transition_summary_df)


###################################################################################



#############Plotting ASR with transition nodes###########
#Load ggtree
library(ggtree)
library(ggtreeExtra,
        lib.loc = "/drives/4tb/Juliet/miniconda3/envs/my_r/lib/R/library/")

p <- ggtree(ASR_models[[1]]$phy) # + geom_tiplab(size = 0.2)

#Change WRP number column name
sxData_rename <- sxData %>% dplyr::rename(`Category of worker reproductive potential` = WRP)
sxData_rename$`Category of worker reproductive potential` <- as.factor(sxData_rename$`Category of worker reproductive potential`)

# Add node labels
# 389 tips, 388 nodes
# tips 1 - 389, nodes 390 onwards
# code says start at 390 and create a sequence of numbers for an additional 387 nodes
ASR_models[[1]]$phy$node.label <- c(seq(390,(390+387))) ## double check


#################################################### 
# SIMPLE TREE WITH LABELLED NODES FOR ADDING SUBFAMILY LABELS
simple_tree <- ggtree(ASR_models[[1]]$phy, 
                     layout = "fan", 
                     branch.length = "none", 
                     open.angle = 25, 
                     size = 0.2) + 
  geom_rootedge(TRUE, size = 0.1) + 
  geom_tiplab(size = 1, offset = 2) +  # Increase tip label size
  geom_nodelab(size = 1)             # Increase node label size

simple_tree


####################################################   


# FIGURE TREE

sxData_rename$`Category of worker reproductive potential` <- factor(sxData_rename$`Category of worker reproductive potential`,
                                                                    levels = c(1, 2, 3, 4), 
                                                                    labels = c("full", "reduced", "males only", "no offspring"))

#basic tree
TEST <- ggtree(ASR_models[[1]]$phy, 
                    layout = "fan", 
                    branch.length = "none", 
                    open.angle = 25, 
                    size = 0.2) + 
  geom_rootedge(TRUE, size = 0.1) #+ 
  #geom_tiplab(size = 1, offset = 2) 

TEST

nicer_tree_rotated <- rotate_tree(TEST, angle=90)

tree_styled <- nicer_tree_rotated +
  geom_fruit(data=sxData_rename, # Data
             geom=geom_tile, # Plots 'tiles' (squares) onto each tip
             width=1.6, # Alters the length of the tiles
             position=position_identityx(hexpand=27), # Adjusts how far away the tiles are from the tree
             mapping=aes(y=animal, fill=`Category of worker reproductive potential`), # Analogous to ggplot aes()
             color = "white", # Colour of outline round the tiles (white makes it look like a gap)
             lwd = 0.2, # Width of line between tiles
             linetype = 1, # Default - other numbers make the line dashes, dotted etc.
             axis.params=list( # Add label to the geom_tile - by adding an x-axis
               axis="x",
               text = " ", # Label to plot
               text.size = 3.5, # Size of text
               hjust = 0, # Adjust position of text relative to the geom_tile
               vjust = 0.8,
               text.angle=360,
               line.colour="white"), # Set to white so axis line is not visible
             offset = 4, # Fine-scale adjustment of space between tree & layers
             pwidth=100) + # Width of whole plot
  scale_fill_manual(values=c("#D33682","#F4A6D3", "#FAD2E1","#FDEDEC")) + # swapped order 17/01/2025
  theme_tree(legend.position = "none") # This line removes the legend

tree_styled

############### ############### ############### ############### ############### ############### ############### 

############### NODE LABELS LOSSES AND GAINS #####################

# Subset the first dataframe from the TransDat_neg_rep_list list of dataframes for rows with particular transition types

table(TransDat_neg_rep_list[[1]]$trans) #model 1 transitions
#WRP_1_to_WRP_1      WRP_1_to_WRP_2    WRP_1_to_WRP_3     WRP_1_to_WRP_4    WRP_2_to_WRP_1    WRP_2_to_WRP_2    WRP_2_to_WRP_3    WRP_3_to_WRP_2 
#.   68                    1                 3                 10                1                      10                   1                   2 
#WRP_3_to_WRP_3.    WRP_4_to_WRP_1     WRP_4_to_WRP_2     WRP_4_to_WRP_3    WRP_4_to_WRP_4 
#.    370              1                        1             24                 284 

# Transitions for losses of worker reproductive potential
one_to_two_loss <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_1_to_WRP_2", trans))
one_to_three_loss <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_1_to_WRP_3", trans))
one_to_four_loss <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_1_to_WRP_4", trans))
two_to_three_loss <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_2_to_WRP_3", trans))
two_to_four_loss <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_2_to_WRP_4", trans)) #empty
three_to_four_loss <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_3_to_WRP_4", trans)) #empty

# Transitions for gains of worker reproductive potential
two_to_one_gain <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_2_to_WRP_1", trans))
three_to_one_gain <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_3_to_WRP_1", trans)) #empty
three_to_two_gain <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_3_to_WRP_2", trans))
four_to_one_gain <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_4_to_WRP_1", trans))
four_to_two_gain <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_4_to_WRP_2", trans))
four_to_three_gain <- subset(TransDat_neg_rep_list[[1]], grepl("WRP_4_to_WRP_3", trans))

# Identify transitional nodes where reproductive capacity decreases (losses)
one_to_two_loss_trans <- one_to_two_loss$V2
one_to_three_loss_trans <- one_to_three_loss$V2
one_to_four_loss_trans <- one_to_four_loss$V2
two_to_three_loss_trans <- two_to_three_loss$V2
two_to_four_loss_trans <- two_to_four_loss$V2
three_to_four_loss_trans <- three_to_four_loss$V2

# Identify transitional nodes where reproductive capacity increases (gains)
two_to_one_gain_trans <- two_to_one_gain$V2
three_to_one_gain_trans <- three_to_one_gain$V2
three_to_two_gain_trans <- three_to_two_gain$V2
four_to_one_gain_trans <- four_to_one_gain$V2
four_to_two_gain_trans <- four_to_two_gain$V2
four_to_three_gain_trans <- four_to_three_gain$V2

# Add transitions for losses of reproductive potential


library(ggplot2)
library(ggtree)

# Add transitions for losses and gains, and include a legend for nodes
tree_trans <- tree_styled + 
  # LOSSES (Mapping color and shape for legend)
  ggtree::geom_point2(aes(subset = node %in% one_to_two_loss_trans,
                          shape = 'Losses of reproductive potential', colour = 'Losses of reproductive potential'),
                      size = 2, fill = 'black', alpha = 0.65) +
  ggtree::geom_point2(aes(subset = node %in% one_to_three_loss_trans,
                          shape = 'Losses of reproductive potential', colour = 'Losses of reproductive potential'),
                      size = 2, fill = 'black', alpha = 0.65) +
  ggtree::geom_point2(aes(subset = node %in% one_to_four_loss_trans,
                          shape = 'Losses of reproductive potential', colour = 'Losses of reproductive potential'),
                      size = 2, fill = 'black', alpha = 0.6) +
  ggtree::geom_point2(aes(subset = node %in% two_to_three_loss_trans,
                          shape = 'Losses of reproductive potential', colour = 'Losses of reproductive potential'),
                      size = 2, fill = 'black', alpha = 0.65) +
  ggtree::geom_point2(aes(subset = node %in% two_to_four_loss_trans,
                          shape = 'Losses of reproductive potential', colour = 'Losses of reproductive potential'),
                      size = 2, fill = 'black', alpha = 0.6) +
  ggtree::geom_point2(aes(subset = node %in% three_to_four_loss_trans,
                          shape = 'Losses of reproductive potential', colour = 'Losses of reproductive potential'),
                      size = 2, fill = 'black', alpha = 0.6) +
  
  # GAINS (Mapping color and shape for legend)
  ggtree::geom_point2(aes(subset = node %in% two_to_one_gain_trans,
                          shape = 'Gains of reproductive potential', colour = 'Gains of reproductive potential'),
                      size = 2, fill = '#D33682', alpha = 0.65) +
  ggtree::geom_point2(aes(subset = node %in% three_to_one_gain_trans,
                          shape = 'Gains of reproductive potential', colour = 'Gains of reproductive potential'),
                      size = 2, fill = '#D33682', alpha = 0.65) +
  ggtree::geom_point2(aes(subset = node %in% three_to_two_gain_trans,
                          shape = 'Gains of reproductive potential', colour = 'Gains of reproductive potential'),
                      size = 2, fill = '#D33682', alpha = 0.65) +
  ggtree::geom_point2(aes(subset = node %in% four_to_one_gain_trans,
                          shape = 'Gains of reproductive potential', colour = 'Gains of reproductive potential'),
                      size = 2, fill = '#D33682', alpha = 0.6) +
  ggtree::geom_point2(aes(subset = node %in% four_to_two_gain_trans,
                          shape = 'Gains of reproductive potential', colour = 'Gains of reproductive potential'),
                      size = 2, fill = '#D33682', alpha = 0.6) +
  ggtree::geom_point2(aes(subset = node %in% four_to_three_gain_trans,
                          shape = 'Gains of reproductive potential', colour = 'Gains of reproductive potential'),
                      size = 2, fill = '#D33682', alpha = 0.6) +
  
  # Fill scale for the tip tiles (reproductive potential)
  scale_fill_manual(values=c("#D33682","#F4A6D3", "#FAD2E1","#FDEDEC"), #swapped order 17/01/2025
                    name = "Category of Worker\nReproductive Potential") +
  
  # Color and shape scale for the node points (transitions)
  scale_colour_manual(name = "Transitions",
                      values = c('Losses of reproductive potential' = 'black', 'Gains of reproductive potential' = '#D33682')) +
  scale_shape_manual(name = "Transitions",
                     values = c('Losses of reproductive potential' = 20, 'Gains of reproductive potential' = 17)) +
  
  # Adding the legend for the plot
  theme(legend.position = "right") # Adjust legend position

tree_trans




################################################################################################
#### ADD SUBFAMILY LABELS



# specify node at base of subfamily clade
tree_clad_lab <- tree_trans + geom_cladelab(node=512, label="Pheidole\n(big-headed ants)", align=TRUE, fontsize = 3, angle="auto",
                                            offset = 2, offset.text = 0.5, textcolor='black', barcolor='black') 

# double check Attine ants
tree_clad_lab <- tree_clad_lab + geom_cladelab(node=563, label="Fungus-growing \nants (including \nAtta and\nAcromyrmex)", align=TRUE, fontsize = 5, angle="auto",
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=453, label="Carebara", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=633, label="Camponotus", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=681, label="Cataglyphis", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=677, label="Formica", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=685, label="Lasius", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')


tree_clad_lab <- tree_clad_lab + geom_cladelab(node=477, label="Temnothorax", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=405, label="Crematogaster", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=711, label="Ponerinae", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=688, label="Dolichoderinae", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=584, label="Solenopsidini \n(including fire ants)", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=602, label="Stenammini", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=624, label="Ectatommini", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=683, label="Oecophylla", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=464, label="Tetramorium", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=462, label="Metapone", align=TRUE, fontsize = 3, angle='auto',
 offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=579, label="Cyphomyrmex", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=701, label="Myrmecia", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=448, label="Cardiocondyla", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=684, label="Plagiolepis", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=583, label="Apterostigma", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')



tree_clad_lab

with_species <- tree_clad_lab + geom_tiplab(size = 1, offset = 2) 
with_species
