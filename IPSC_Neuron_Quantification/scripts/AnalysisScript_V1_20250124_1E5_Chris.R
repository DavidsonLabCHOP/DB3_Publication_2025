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


measurements <- read_delim("/Users/christopherfluta/Documents/DB3_IPSC_Neurons_20250110/Measurements/measurements_New.tsv", 
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
  
  #  For each "Cell" in the dataframe, determine the mean intensity level of blue
  Blue_Level <- workingImage$`Nucleus: Blue mean`
  
  # For each "Cell" in the dataframe, determine the nucleus size
  Nuc_Size <- workingImage$`Nucleus: Area`
  
  # Create a dataframe with results
  perCellPunctaCountsTable <- data.frame(Cell_IDs, Red_Level, Green_Level, Blue_Level, Nuc_Size)
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
      x$Nuc_Size >= 700 &
      x$Blue_Level < 40 &
      (x$Blue_Level / x$Nuc_Size) < 0.05,
  ]
  
  # List settings
  # Size parameter has a good effect on excluding extra cells 
  # PPT with example images fromm 300/400 size vs 600 size
  # 400 vs 600 plot side by side with same y scale
  
  DB3_Positive_Cells <- x[x$Green_Level >= 10,]
  MAP2_DB3_Positive_Cells <- intersect(MAP2_Positive_Cells$Cell_IDs, DB3_Positive_Cells$Cell_IDs)
  MAP2_Positive_Count <- length(MAP2_Positive_Cells$Cell_IDs)
  DB3_Positive_Count <- length(DB3_Positive_Cells$Cell_IDs)
  MAP2_DB3_Positive_Cells <- length(MAP2_DB3_Positive_Cells)
  
  # Return a named vector for clarity
  return(c(MAP2_Positive_Count = MAP2_Positive_Count,
           DB3_Positive_Count = DB3_Positive_Count,
           MAP2_DB3_Positive_Cells = MAP2_DB3_Positive_Cells))
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

#write.csv(tallys_for_plot, file = "/Users/christopherfluta/Documents/DB3_IPSC_Neurons_20250110/Figures/BarplotTable.csv")

melted_tallys_for_plot <- melt(tallys_for_plot)

# Graph
ggplot(melted_tallys_for_plot, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("IPSC Treatment Condition") + 
  ylab("Number of Positive Cells")

melted_tallys_for_plot$Treatment <- ifelse(
  grepl("^AAV1", melted_tallys_for_plot$Sample), "AAV1",
  ifelse(grepl("^DB3", melted_tallys_for_plot$Sample), "DB3", "Other")
)

melted_tallys_for_plot$Dosage <- sub(".*_(1E[0-9])", "\\1", melted_tallys_for_plot$Sample)

melted_tallys_for_plot$Sample <- factor(
  melted_tallys_for_plot$Sample,
  levels = c("AAV1_1E3", "DB3_1E3", "AAV1_1E4", "DB3_1E4", "AAV1_1E5", "DB3_1E5")
)

colnames(melted_tallys_for_plot) <- c("Sample", "Variable", "value","Treatment","Dosage")

dosage_labels <- data.frame(
  x = c(1.5, 3.5, 5.5), # Positions for centered labels
  label = c("1E+03", "1E+04", "1E+05") # Labels for each dosage group
)

library(ggplot2)

melted_tallys_for_plot$Sample <- factor(
  melted_tallys_for_plot$Sample,
  levels = c(
    "AAV1_1E3", "DB3_1E3", "Spacer1", # Add a spacer after 1E3
    "AAV1_1E4", "DB3_1E4", "Spacer2", # Add a spacer after 1E4
    "AAV1_1E5", "DB3_1E5"             # No spacer after the last group
  )
)

# Add dummy rows for spacers with NA values
dummy_rows <- data.frame(
  Sample = c("Spacer1", "Spacer2"),
  Variable = NA,
  value = NA,
  Treatment = NA,
  Dosage = NA
)
melted_tallys_for_plot <- rbind(melted_tallys_for_plot, dummy_rows)

# Plot with spacers
ggplot(melted_tallys_for_plot, aes(x = Sample, y = value, fill = Variable)) +
  # Add bars for the mean
  stat_summary(fun = "mean", geom = "bar", position = "dodge", alpha = 0.7, aes(color = Treatment)) +
  # Add error bars for SEM
  stat_summary(
    fun.data = function(x) {
      data.frame(
        y = mean(x),
        ymin = mean(x) - sd(x) / sqrt(length(x)),
        ymax = mean(x) + sd(x) / sqrt(length(x))
      )
    },
    geom = "errorbar",
    position = position_dodge(0.9),
    width = 0.25
  ) +
  ylim(-5, 125) +
  # Custom fill colors for bars
  scale_fill_manual(
    values = c("lightgrey", "grey50", "grey0")
  ) +
  # Custom color scale for outlines
  scale_color_manual(
    values = c("AAV1" = "black", "DB3" = "white")
  ) +
  # Adjust x-axis
  scale_x_discrete(
    labels = c(
      "AAV1_1E3" = "1E+03", "DB3_1E3" = "", "Spacer1" = "", # Labels for 1E3 group
      "AAV1_1E4" = "1E+04", "DB3_1E4" = "", "Spacer2" = "", # Labels for 1E4 group
      "AAV1_1E5" = "1E+05", "DB3_1E5" = ""                 # Labels for 1E5 group
    )
  ) +
  # Add labels and theme
  xlab("Dosage") +
  ylab("Mean Positive Cells Per Image") +
  theme_minimal() +
  theme_classic() +  # Use classic theme to remove grey background and gridlines
  theme(
    axis.line = element_line(color = "black"), # Ensure solid black axis lines
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.text.x = element_text(angle = 0, size = 9, hjust = -.4, vjust = 11), # Add angled x-axis text
    legend.position = "right" # Adjust legend position if needed
  )

ggsave(
  filename = "nuc700_300dpi_y125.pdf",        # Specify the output file name
  plot = last_plot(),           # Specify the plot to save (or replace with your plot object)
  device = "pdf",               # Set the file format to PDF
  dpi = 300,                    # Set the resolution to 300 DPI
  width = 8,                    # Set the width of the plot in inches
  height = 6,                   # Set the height of the plot in inches
  units = "in"                  # Specify the unit of measurement for width and height
)