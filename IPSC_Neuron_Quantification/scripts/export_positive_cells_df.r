library(readr)
library(vioplot)
library(ggplot2)
library(dplyr)
#library(tidyr)
#library(reshape2)

################
# Set Thresholds
# Darp32 positive puncta threshold 
#D32Threshold <- 20 
#mTFP1Threshold <- 20
#
################


measurements <- read_delim("/Users/christopherfluta/Documents/DB3_IPSC_Neurons_20250110/Measurements/measurements_Chris.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

Image_Names <- levels(as.factor(measurements$Image))

# Split datasets into 1E3 and 1E4 categories
ImageNames_1E3 <- Image_Names[grep("_1E3 ", Image_Names)]
ImageNames_1E4 <- Image_Names[grep("_1E4 ", Image_Names)]
ImageNames_1E5 <- Image_Names[grep("_1E5 ", Image_Names)]

# Separate out by sample and dose
AAV1_1E3 <- ImageNames_1E3[grep("AAV1-mRuby3g", ImageNames_1E3)]
AAV1_1E4 <- ImageNames_1E4[grep("AAV1-mRuby3g", ImageNames_1E4)]
AAV1_1E5 <- ImageNames_1E5[grep("AAV1-mRuby3g", ImageNames_1E5)]

DB3_1E3 <- ImageNames_1E3[grep("DB3-mTFP1g", ImageNames_1E3)]
DB3_1E4 <- ImageNames_1E4[grep("DB3-mTFP1g", ImageNames_1E4)]
DB3_1E5 <- ImageNames_1E5[grep("DB3-mTFP1g", ImageNames_1E5)]

SampleList <- list(AAV1_1E3,AAV1_1E4,AAV1_1E5,
                   DB3_1E3,DB3_1E4,DB3_1E5)

#For testing
#imageID <- "240923 MSH3 DARPP32 in mouse.lif - 270924 WT control 1 field 1 40x 520 MSH3 620 DARP32_Processed001.tif" 

# function to split output data into a list of dataframes, one for each image.
processImage_fun <- function(imageID){
  workingImage <- measurements[grep(imageID, measurements$Image),]
  
  return(workingImage)
}

# Run function to parse images into a list
AAV1_1E3_images_list <- lapply(AAV1_1E3, processImage_fun)
AAV1_1E4_images_list <- lapply(AAV1_1E4, processImage_fun)
AAV1_1E5_images_list <- lapply(AAV1_1E5, processImage_fun)
DB3_1E3_images_list <- lapply(DB3_1E3, processImage_fun)
DB3_1E4_images_list <- lapply(DB3_1E4, processImage_fun)
DB3_1E5_images_list <- lapply(DB3_1E5, processImage_fun)

# Create a function to collect cell abundance information
cell_count_fun <- function(image_dataframe){
  workingImage <- as.data.frame(image_dataframe)
  cells_df <- workingImage[workingImage$`Object type` == "Cell", ]
  #cells_num <- dim(cells_df)
  cells_num <- length(cells_df$`Object type`)
  return(cells_num)
}

# Get cell abundance from each image 
AAV1_1E3_cell_count <- unlist(lapply(AAV1_1E3_images_list, cell_count_fun))
AAV1_1E4_cell_count <- unlist(lapply(AAV1_1E4_images_list, cell_count_fun))
AAV1_1E5_cell_count <- unlist(lapply(AAV1_1E5_images_list, cell_count_fun))
DB3_1E3_cell_count <- unlist(lapply(DB3_1E3_images_list, cell_count_fun))
DB3_1E4_cell_count <- unlist(lapply(DB3_1E4_images_list, cell_count_fun))
DB3_1E5_cell_count <- unlist(lapply(DB3_1E5_images_list, cell_count_fun))


# Create a function to assign cell IDs to cells and detections
assign_cell_ids <- function(image_dataframe) {
  # Read the image dataframe
  data <- image_dataframe
  
  # Initialize Cell_ID counter
  cell_id_counter <- 0
  
  # Create a vector to store Cell_IDs
  cell_ids <- integer(nrow(data))
  
  # Iterate over each row in the data frame
  for (i in seq_len(nrow(data))) {
    if (data$`Object type`[i] == "Cell") {
      # Increment Cell_ID counter when a Cell object is encountered
      cell_id_counter <- cell_id_counter + 1
    }
    # Assign current Cell_ID to both Cell and Detection objects
    cell_ids[i] <- cell_id_counter
  }
  
  # Add the Cell_ID column to the data frame
  data$Name <- cell_ids
  
  return(data)
}

AAV1_1E3_images_list_cellIDs <- lapply(AAV1_1E3_images_list, assign_cell_ids)
AAV1_1E4_images_list_cellIDs <- lapply(AAV1_1E4_images_list, assign_cell_ids)
AAV1_1E5_images_list_cellIDs <- lapply(AAV1_1E5_images_list, assign_cell_ids)
DB3_1E3_images_list_cellIDs <- lapply(DB3_1E3_images_list, assign_cell_ids)
DB3_1E4_images_list_cellIDs <- lapply(DB3_1E4_images_list, assign_cell_ids)
DB3_1E5_images_list_cellIDs <- lapply(DB3_1E5_images_list, assign_cell_ids)

# Create a function to count the number of transduced neurons in each image
punctaCount_fun <- function(image_dataframe){
  
  #workingImage <- as.data.frame(miMSH311_115G_3_images_list_cellIDs)
  workingImage <- as.data.frame(image_dataframe)
  
  # split dataframe into a list of dataframes one for each cell Index the list by cell ID and include only detections in each table.
  # Get vector of cell IDs
  Cell_IDs <- workingImage$Name
  
  # Subset dataframe to include only Cells
  #detections_df <- workingImage[workingImage$`Object.type` == "Cell", ]
  
  # Depricated as we are using cells here instead of detections
  # # Split detections_df into a list on Cell_IDs
  # detections_perCell_List <- lapply(Cell_IDs, function(x){
  #   detection_perCell_df <- detections_df[detections_df$Name == x,]
  #   return(detection_perCell_df)
  # })
  # 
  # # Add Cell_ID names to the list
  # names(detections_perCell_List) <- Cell_IDs
  # 
  
  
  # For each "Cell" in the dataframe, determine the mean intensity level of red
  Red_Level <- workingImage$`Nucleus: Red mean`
  
  # For each "Cell" in the dataframe, determine the mean intensity level of green
  Green_Level <- workingImage$`Nucleus: Green mean`
  
  # For each "Cell" in the dataframe, determine the nucleus size
  Nuc_Size <- workingImage$`Nucleus: Area`
  
  Blue_Level <- workingImage$`Nucleus: Blue mean`      # Add Blue Level
  
  Object_ID <- workingImage$`Object ID`
  
  # Create a dataframe with results
  perCellPunctaCountsTable <- data.frame(Cell_IDs, Red_Level, Green_Level, Blue_Level, Nuc_Size,Object_ID)
  return(perCellPunctaCountsTable)
}

# Test the above function on one image.
#punctaCount_fun(AAV1_1E3_images_list_cellIDs[[1]])

# Run the puncata counting function on all images in list
AAV1_1E3_puncatCount_perImage <- lapply(AAV1_1E3_images_list_cellIDs, punctaCount_fun)
AAV1_1E4_puncatCount_perImage <- lapply(AAV1_1E4_images_list_cellIDs, punctaCount_fun)
AAV1_1E5_puncatCount_perImage <- lapply(AAV1_1E5_images_list_cellIDs, punctaCount_fun)
DB3_1E3_puncatCount_perImage <- lapply(DB3_1E3_images_list_cellIDs, punctaCount_fun)
DB3_1E4_puncatCount_perImage <- lapply(DB3_1E4_images_list_cellIDs, punctaCount_fun)
DB3_1E5_puncatCount_perImage <- lapply(DB3_1E5_images_list_cellIDs, punctaCount_fun)

# Visualize the distributions
par(mfrow=c(4, 2)) # use a 3x2 (rows x columns) layout

lapply(DB3_1E4_puncatCount_perImage, vioplot) # call plot for each list element

par(mfrow=c(1, 1)) # reset layout

# Create a function to calculate the number of red positive cells, the number of green positive cells, and the number of red green double positive cells
tallyCells_Fun <- function(x){
  ### Interesect method
  #MAP2_Positive_Cells_1 <- x[x$Red_Level >= 10,]
  #MAP2_Positive_Cells_2 <- x[x$Nuc_Size >= 600,]
  #MAP2_Positive_Cells_3 <- x[x$Blue_Level < 20]
  #MAP2_Positive_Cells_4 <- x[(x$Blue_Level / x$Nuc_Size) < 0.5, ] 
  
  #MAP2_Positive_Cells <- intersect(MAP2_Positive_Cells_1,MAP2_Positive_Cells_2,MAP2_Positive_Cells_3,MAP2_Positive_Cells_4)
  
  # One operation method
  MAP2_Positive_Cells <- x[
    x$Red_Level >= 7.5 &
      x$Nuc_Size >= 800 &
      x$Blue_Level < 40 &
      (x$Blue_Level / x$Nuc_Size) < 0.05,
  ]
  
  DB3_Positive_Cells <- x[x$Green_Level >= 10,]
  MAP2_DB3_Positive_Cells <- intersect(MAP2_Positive_Cells$Cell_IDs, DB3_Positive_Cells$Cell_IDs)
  MAP2_Positive_Count <- length(MAP2_Positive_Cells$Cell_IDs)
  DB3_Positive_Count <- length(DB3_Positive_Cells$Cell_IDs)
  MAP2_DB3_Positive_Cells <- length(MAP2_DB3_Positive_Cells)
  
  # Return a named vector for clarity
  return(MAP2_Positive_Cells)
}

# Use lapply to apply the function and do.call with rbind to create a data frame
AAV1_1E3_CellCounts <- do.call(rbind, lapply(AAV1_1E3_puncatCount_perImage, tallyCells_Fun))
# Convert to a data frame
AAV1_1E3_CellCounts <- as.data.frame(AAV1_1E3_CellCounts)
AAV1_1E3_CellCounts["Sample"] <- "AAV1_1E3"


# Use lapply to apply the function and do.call with rbind to create a data frame
AAV1_1E4_CellCounts <- do.call(rbind, lapply(AAV1_1E4_puncatCount_perImage, tallyCells_Fun))
# Convert to a data frame
AAV1_1E4_CellCounts <- as.data.frame(AAV1_1E4_CellCounts)
AAV1_1E4_CellCounts["Sample"] <- "AAV1_1E4"

# Use lapply to apply the function and do.call with rbind to create a data frame
AAV1_1E5_CellCounts <- do.call(rbind, lapply(AAV1_1E5_puncatCount_perImage, tallyCells_Fun))
# Convert to a data frame
AAV1_1E5_CellCounts <- as.data.frame(AAV1_1E5_CellCounts)
AAV1_1E5_CellCounts["Sample"] <- "AAV1_1E5"

# Use lapply to apply the function and do.call with rbind to create a data frame
DB3_1E3_CellCounts <- do.call(rbind, lapply(DB3_1E3_puncatCount_perImage, tallyCells_Fun))
# Convert to a data frame
DB3_1E3_CellCounts <- as.data.frame(DB3_1E3_CellCounts)
DB3_1E3_CellCounts["Sample"] <- "DB3_1E3"

# Use lapply to apply the function and do.call with rbind to create a data frame
DB3_1E4_CellCounts <- do.call(rbind, lapply(DB3_1E4_puncatCount_perImage, tallyCells_Fun))
# Convert to a data frame
DB3_1E4_CellCounts <- as.data.frame(DB3_1E4_CellCounts)
DB3_1E4_CellCounts["Sample"] <- "DB3_1E4"

# Use lapply to apply the function and do.call with rbind to create a data frame
DB3_1E5_CellCounts <- do.call(rbind, lapply(DB3_1E5_puncatCount_perImage, tallyCells_Fun))
# Convert to a data frame
DB3_1E5_CellCounts <- as.data.frame(DB3_1E5_CellCounts)
DB3_1E5_CellCounts["Sample"] <- "DB3_1E5"

library(reshape2)
tallys_for_plot <- rbind(AAV1_1E3_CellCounts,AAV1_1E4_CellCounts,AAV1_1E5_CellCounts,
                         DB3_1E3_CellCounts,DB3_1E4_CellCounts,DB3_1E5_CellCounts)

#tallys_for_plot <- tallys_for_plot[6]
write.csv(tallys_for_plot, file = "/Users/christopherfluta/Documents/DB3_IPSC_Neurons_20250110/Figures/MAP2_Pos_cells_full_800.csv",
          row.names = FALSE)
