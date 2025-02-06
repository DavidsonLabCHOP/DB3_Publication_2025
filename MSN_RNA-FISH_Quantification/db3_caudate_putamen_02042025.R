# Load necessary libraries
library(dplyr)
library(writexl)

# Define thresholds
#ch4_threshold_c <- 2
#ch4_threshold_a_b <- 10
ch2_drd_puncta_threshold <- 45
ch3_drd_puncta_threshold <- 30
ch5_neun_threshold <- 10

# Process spreadsheet function
process_spreadsheet <- function(filepath, ch4_threshold) {
  # Load the data
  data <- read.table(filepath, header = TRUE, sep = "\t")
  
  # Filter rows where Object.type is "Cell"
  data_filtered <- data[data$Object.type == "Cell", ]
  
  # Extract animal ID and region
  data_filtered$Animal_ID <- sapply(strsplit(data_filtered$Image, " "), `[`, 2)
  data_filtered$Region <- ifelse(grepl("caudate", data_filtered$Image, ignore.case = TRUE), "Caudate",
                                 ifelse(grepl("putamen", data_filtered$Image, ignore.case = TRUE), "Putamen", NA))
  
  # Filter out rows without a defined region
  data_filtered <- data_filtered[!is.na(data_filtered$Region), ]
  
  # Check if ch5 exists in the dataset
  ch5_present <- "Subcellular..Channel.5..Num.spots.estimated" %in% colnames(data_filtered)
  
  # Initialize an empty data frame for results
  results <- data.frame(
    Image = character(),
    Animal_ID = character(),
    Region = character(),
    Total_Cells = integer(),
    Cells_Transduced = integer(),
    Total_Non_MSNS = integer(),
    Non_MSNS_Transduced = integer(),
    Total_MSNs = integer(),
    MSNs_Transduced = integer(),
    Channel2_Cells = integer(),
    Transduced_Channel2 = integer(),
    Channel3_Cells = integer(),
    Transduced_Channel3 = integer(),
    Channel2_and_3_Cells = integer(),
    Transduced_Channel2_and_3 = integer(),
    Total_Neun = integer(),
    Neun_Transduced = integer(),
    stringsAsFactors = FALSE
  )
  
  # Process each image
  unique_images <- unique(data_filtered$Image)
  for (image in unique_images) {
    image_data <- data_filtered[data_filtered$Image == image, ]
    animal_id <- unique(image_data$Animal_ID)
    region <- unique(image_data$Region)
    
    total_cells <- 0
    cells_transduced <- 0
    total_non_msns <- 0
    non_msns_transduced <- 0
    total_msns <- 0
    msns_transduced <- 0
    channel2_cells <- 0
    channel3_cells <- 0
    channel2_and_3_cells <- 0
    transduced_channel2 <- 0
    transduced_channel3 <- 0
    transduced_channel2_and_3 <- 0
    total_neun <- 0
    neun_transduced <- 0
    
    for (i in 1:nrow(image_data)) {
      total_cells <- total_cells + 1
      
      channel2_value <- as.numeric(image_data[i, "Subcellular..Channel.2..Num.spots.estimated"])
      channel3_value <- as.numeric(image_data[i, "Subcellular..Channel.3..Num.spots.estimated"])
      channel4_value <- as.numeric(image_data[i, "Subcellular..Channel.4..Num.spots.estimated"])
      channel5_value <- if (ch5_present) as.numeric(image_data[i, "Subcellular..Channel.5..Num.spots.estimated"]) else NA
      
      # Transduced check
      if (!is.na(channel4_value) && channel4_value >= ch4_threshold) {
        cells_transduced <- cells_transduced + 1
      }
      
      # Neun classification
      if (ch5_present && !is.na(channel5_value) && channel5_value >= ch5_neun_threshold) {
        total_neun <- total_neun + 1
        if (!is.na(channel4_value) && channel4_value >= ch4_threshold) {
          neun_transduced <- neun_transduced + 1
        }
      }
      
      # Default to non-MSN, override if MSN conditions are met
      is_msn <- FALSE
      
      # Check if cell meets DRD1 or DRD2 threshold
      if (!is.na(channel2_value) && !is.na(channel3_value)) {
        total_puncta <- channel2_value + channel3_value
        channel2_pct <- channel2_value / total_puncta
        channel3_pct <- channel3_value / total_puncta
        
        if (channel2_value >= ch2_drd_puncta_threshold && channel2_pct > 0.6) {
          # DRD1 cell
          channel2_cells <- channel2_cells + 1
          total_msns <- total_msns + 1
          is_msn <- TRUE
          if (channel4_value >= ch4_threshold) {
            transduced_channel2 <- transduced_channel2 + 1
            msns_transduced <- msns_transduced + 1
          }
        } else if (channel3_value >= ch3_drd_puncta_threshold && channel3_pct > 0.6) {
          # DRD2 cell
          channel3_cells <- channel3_cells + 1
          total_msns <- total_msns + 1
          is_msn <- TRUE
          if (channel4_value >= ch4_threshold) {
            transduced_channel3 <- transduced_channel3 + 1
            msns_transduced <- msns_transduced + 1
          }
        } else if (channel2_value >= ch2_drd_puncta_threshold && channel3_value >= ch3_drd_puncta_threshold) {
          # Double positive MSN
          channel2_and_3_cells <- channel2_and_3_cells + 1
          total_msns <- total_msns + 1
          is_msn <- TRUE
          if (channel4_value >= ch4_threshold) {
            transduced_channel2_and_3 <- transduced_channel2_and_3 + 1
            msns_transduced <- msns_transduced + 1
          }
        }
      }
      
      # If not classified as MSN, classify as Non-MSN
      if (!is_msn) {
        total_non_msns <- total_non_msns + 1
        if (!is.na(channel4_value) && channel4_value >= ch4_threshold) {
          non_msns_transduced <- non_msns_transduced + 1
        }
      }
    }
    
    # Append results
    results <- rbind(results, data.frame(
      Image = image,
      Animal_ID = animal_id,
      Region = region,
      Total_Cells = total_cells,
      Cells_Transduced = cells_transduced,
      Total_Non_MSNS = total_non_msns,
      Non_MSNS_Transduced = non_msns_transduced,
      Total_MSNs = total_msns,
      MSNs_Transduced = msns_transduced,
      Channel2_Cells = channel2_cells,
      Transduced_Channel2 = transduced_channel2,
      Channel3_Cells = channel3_cells,
      Transduced_Channel3 = transduced_channel3,
      Channel2_and_3_Cells = channel2_and_3_cells,
      Transduced_Channel2_and_3 = transduced_channel2_and_3,
      Total_Neun = if (ch5_present) total_neun else NA,
      Neun_Transduced = if (ch5_present) neun_transduced else NA,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

# Process spreadsheets
results_a <- process_spreadsheet("/Users/christopherfluta/Documents/Neun_drd_measurements_14.tsv", 10)
results_b <- process_spreadsheet("/Users/christopherfluta/Downloads/custom_measurements_14all.tsv", 10)
results_c <- process_spreadsheet("/Volumes/T7/db3injected-1nhp-rh16070521-mtfp1-drd1-drd2-230317/measurements_DRD50_Ch4_83.tsv", 10)

combined_results <- bind_rows(results_a, results_b, results_c)

print(combined_results)


# Summarize and calculate percentages
percentages_summary <- combined_results %>%
  group_by(Animal_ID, Region) %>%
  summarize(
    Total_Cells = sum(Total_Cells),
    Cells_Transduced = sum(Cells_Transduced),
    Total_Non_MSNS = sum(Total_Non_MSNS),
    Non_MSNS_Transduced = sum(Non_MSNS_Transduced),
    Total_MSNs = sum(Total_MSNs),
    MSNs_Transduced = sum(MSNs_Transduced),
    Channel2_Cells = sum(Channel2_Cells),
    Transduced_Channel2 = sum(Transduced_Channel2),
    Channel3_Cells = sum(Channel3_Cells),
    Transduced_Channel3 = sum(Transduced_Channel3),
    Channel2_and_3_Cells = sum(Channel2_and_3_Cells),
    Transduced_Channel2_and_3 = sum(Transduced_Channel2_and_3),
    Total_Neun = sum(Total_Neun, na.rm = TRUE),
    Neun_Transduced = sum(Neun_Transduced, na.rm = TRUE),
    percentage_transduced_all = Cells_Transduced / Total_Cells * 100,
    percentage_transduced_non_msns = Non_MSNS_Transduced / Total_Non_MSNS * 100,
    percentage_transduced_msns = MSNs_Transduced / Total_MSNs * 100,
    percentage_transduced_ch2 = Transduced_Channel2 / Channel2_Cells * 100,
    percentage_transduced_ch3 = Transduced_Channel3 / Channel3_Cells * 100,
    percentage_transduced_ch2_and_3 = Transduced_Channel2_and_3 / Channel2_and_3_Cells * 100,
    percentage_transduced_neun = ifelse(Total_Neun > 0, Neun_Transduced / Total_Neun * 100, NA),
    .groups = "drop"
  )


print(percentages_summary)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(writexl)

# Drop animal 191393
filtered_results <- combined_results %>%
  filter(Animal_ID != "191393")

#filtered_results <- filtered_results %>%
#  filter(!(Animal_ID == "15071281" & Region == "Caudate"))

# Calculate percentages for individual animals
individual_summary <- filtered_results %>%
  group_by(Animal_ID, Region) %>%
  summarize(
    Total_Cells = sum(Total_Cells),
    Cells_Transduced = sum(Cells_Transduced),
    Total_Non_MSNS = sum(Total_Non_MSNS),
    Non_MSNS_Transduced = sum(Non_MSNS_Transduced),
    Total_MSNs = sum(Total_MSNs),
    MSNs_Transduced = sum(MSNs_Transduced),
    Channel2_Cells = sum(Channel2_Cells),
    Transduced_Channel2 = sum(Transduced_Channel2),
    Channel3_Cells = sum(Channel3_Cells),
    Transduced_Channel3 = sum(Transduced_Channel3),
    Channel2_and_3_Cells = sum(Channel2_and_3_Cells),
    Transduced_Channel2_and_3 = sum(Transduced_Channel2_and_3),
    Total_Neun = sum(Total_Neun, na.rm = TRUE),
    Neun_Transduced = sum(Neun_Transduced, na.rm = TRUE),
    percentage_transduced_all = Cells_Transduced / Total_Cells * 100,
    percentage_transduced_non_msns = Non_MSNS_Transduced / Total_Non_MSNS * 100,
    percentage_transduced_msns = MSNs_Transduced / Total_MSNs * 100,
    percentage_transduced_ch2 = Transduced_Channel2 / Channel2_Cells * 100,
    percentage_transduced_ch3 = Transduced_Channel3 / Channel3_Cells * 100,
    percentage_transduced_ch2_and_3 = Transduced_Channel2_and_3 / Channel2_and_3_Cells * 100,
    percentage_transduced_neun = ifelse(Total_Neun > 0, Neun_Transduced / Total_Neun * 100, NA),
    .groups = "drop"
  )


filtered_results <- filtered_results %>%
  filter(!(Animal_ID == "RH16070521"))

# Calculate combined percentages per region for GP animals
joined_summary <- filtered_results %>%
  group_by(Region) %>%
  summarize(
    Total_Cells = sum(Total_Cells),
    Cells_Transduced = sum(Cells_Transduced),
    Total_Non_MSNS = sum(Total_Non_MSNS),
    Non_MSNS_Transduced = sum(Non_MSNS_Transduced),
    Total_MSNs = sum(Total_MSNs),
    MSNs_Transduced = sum(MSNs_Transduced),
    Channel2_Cells = sum(Channel2_Cells),
    Transduced_Channel2 = sum(Transduced_Channel2),
    Channel3_Cells = sum(Channel3_Cells),
    Transduced_Channel3 = sum(Transduced_Channel3),
    Channel2_and_3_Cells = sum(Channel2_and_3_Cells),
    Transduced_Channel2_and_3 = sum(Transduced_Channel2_and_3),
    Total_Neun = sum(Total_Neun, na.rm = TRUE),
    Neun_Transduced = sum(Neun_Transduced, na.rm = TRUE),
    percentage_transduced_all = Cells_Transduced / Total_Cells * 100,
    percentage_transduced_non_msns = Non_MSNS_Transduced / Total_Non_MSNS * 100,
    percentage_transduced_msns = MSNs_Transduced / Total_MSNs * 100,
    percentage_transduced_ch2 = Transduced_Channel2 / Channel2_Cells * 100,
    percentage_transduced_ch3 = Transduced_Channel3 / Channel3_Cells * 100,
    percentage_transduced_ch2_and_3 = Transduced_Channel2_and_3 / Channel2_and_3_Cells * 100,
    percentage_transduced_neun = ifelse(Total_Neun > 0, Neun_Transduced / Total_Neun * 100, NA),
    .groups = "drop"
  )


# Save results
write_xlsx(
  list(
    "Combined Results" = combined_results,
    "Animal Percentages Summary" = percentages_summary,
    "Joined GP Summary" = joined_summary
  ),
  "Transduction_Percentages_Summary.xlsx"
)

# Reshape data for plotting
metric_labels <- c(
  "percentage_transduced_all" = "All Cells",
  "percentage_transduced_non_msns" = "Non-MSN Cells",
  "percentage_transduced_neun" = "Neun Positive Cells",
  "percentage_transduced_msns" = "All MSNs",
  "percentage_transduced_ch2" = "DRD1 Positive MSNs",
  "percentage_transduced_ch3" = "DRD2 Positive MSNs",
  "percentage_transduced_ch2_and_3" = "Double Positive MSNs"
)

# Filter rows with valid values for plotting
joined_summary_long <- joined_summary %>%
  pivot_longer(
    cols = starts_with("percentage"),
    names_to = "Metric",
    values_to = "Percentage"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Percentage = ifelse(is.na(Percentage), 0, Percentage) # Replace NA with 0
  )

individual_summary_long <- individual_summary %>%
  pivot_longer(
    cols = starts_with("percentage"),
    names_to = "Metric",
    values_to = "Percentage"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Percentage = ifelse(is.na(Percentage), 0, Percentage) # Replace NA with 0
  )

individual_summary_long <- individual_summary_long %>%
  filter(!(Animal_ID == "RH16070521" & Metric == "Neun Positive Cells"))

# Define position dodge for bars and points
dodge_position <- position_dodge(width = 0.9)

plot_gp <- ggplot() +
  # Bars: Combined percentages
  geom_bar(
    data = joined_summary_long,
    aes(x = Metric, y = Percentage, fill = Region, group = Region),
    stat = "identity",
    position = dodge_position,
    alpha = 1
  ) +
  # Points: Individual animal percentages
  geom_point(
    data = individual_summary_long,
    aes(x = Metric, y = Percentage, group = Region),
    position = dodge_position,
    size = 2,
    shape = 21,
    color = "black",
    fill = "black",
    show.legend = FALSE
  ) +
  labs(
    title = "Percentage Transduction in GP Animals",
    x = "Cell Type",
    y = "Percent Transduction"
  ) +
  scale_fill_manual(values = c("#6754E2", "#7c7a8c")) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


#print(plot_gp)


### Begin calculations of transduction with proper RH16070521 files ###

RH16070521_Putamen_Results <- readRDS(file = "/Users/christopherfluta/Downloads/RH16070521_Putamen_Results_Table_NEW_TFP50_min2.rds")
RH16070521_Caudate_Results <- readRDS(file = "/Users/christopherfluta/Downloads/RH16070521_Caudate_Results_Table_NEW_TFP50_min2.rds")

# Filter the relevant columns for each dataset
RH16070521_Caudate_Results_filtered <- RH16070521_Caudate_Results[, 1:14]
RH16070521_Putamen_Results_filtered <- RH16070521_Putamen_Results[, 1:14]

# Calculate column sums for each dataset
caudate_sums <- colSums(RH16070521_Caudate_Results_filtered, na.rm = TRUE)
putamen_sums <- colSums(RH16070521_Putamen_Results_filtered, na.rm = TRUE)

cell_counts <- bind_rows(caudate_sums,putamen_sums)

library(tibble)

# Add row identifiers and convert to row names
cell_counts <- bind_rows(caudate_sums, putamen_sums) %>%
  mutate(Dataset = c("Caudate", "Putamen")) %>%
  column_to_rownames("Dataset")

# View the result
print(cell_counts)

library(dplyr)

# Add percentage columns to the DataFrame
cell_counts <- cell_counts %>%
  mutate(
    all_cells_pct_transduced = TealSum / Dapi * 100,
    non_msmns_pct_transduced = DapiTeal / DapiOnlySum * 100,
    msns_pct_transduced = (DapiMagentaTeal+DapiGreenTeal+DapiGreenMagentaTeal) / (MagentaSum+GreenSum+GreenMagentaSum) * 100,
    drd1_pct_transduced = DapiGreenTeal / GreenSum * 100,
    drd2_pct_transduced = DapiMagentaTeal / MagentaSum * 100,
    drd1_drd2_pct_transduced = DapiGreenMagentaTeal / GreenMagentaSum * 100
  )

# View the result
#print(cell_counts)
cell_counts_subset <- cell_counts[,15:20]

colnames(cell_counts_subset) <- c("percentage_transduced_all", "percentage_transduced_non_msns",
                                  "percentage_transduced_msns", "percentage_transduced_ch2", "percentage_transduced_ch3", "percentage_transduced_ch2_and_3")

# Subset the columns to match the relevant columns in individual_summary
columns_to_replace <- colnames(cell_counts_subset)

# Replace the last two rows with values from cell_counts_subset
individual_summary[5:6, columns_to_replace] <- as.data.frame(cell_counts_subset)

cell_counts_for_joined_animals <- cell_counts[,1:14]


joined_summary <- joined_summary %>%
  mutate(
    Total_Cells = if_else(Region == "Caudate", Total_Cells + 849, 
                          if_else(Region == "Putamen", Total_Cells + 1402, Total_Cells)),
    
    Cells_Transduced = if_else(Region == "Caudate", Cells_Transduced + 259, 
                               if_else(Region == "Putamen", Cells_Transduced + 368, Cells_Transduced)),
    
    Total_Non_MSNS = if_else(Region == "Caudate", Total_Non_MSNS + 382, 
                             if_else(Region == "Putamen", Total_Non_MSNS + 662, Total_Non_MSNS)),
    
    Non_MSNS_Transduced = if_else(Region == "Caudate", Non_MSNS_Transduced + 47, 
                                  if_else(Region == "Putamen", Non_MSNS_Transduced + 92, Non_MSNS_Transduced)),
    
    Total_MSNs = if_else(Region == "Caudate", Total_MSNs + 198+220+49, 
                         if_else(Region == "Putamen", Total_MSNs + 325+295+120, Total_MSNs)),
    
    MSNs_Transduced = if_else(Region == "Caudate", MSNs_Transduced + 84+114+14, 
                              if_else(Region == "Putamen", MSNs_Transduced + 79+177+20, MSNs_Transduced)),
    
    Channel2_Cells = if_else(Region == "Caudate", Channel2_Cells + 198, 
                             if_else(Region == "Putamen", Channel2_Cells + 325, Channel2_Cells)),
    
    Transduced_Channel2 = if_else(Region == "Caudate", Transduced_Channel2 + 114, 
                                  if_else(Region == "Putamen", Transduced_Channel2 + 177, Transduced_Channel2)),
    
    Channel3_Cells = if_else(Region == "Caudate", Channel3_Cells + 220, 
                             if_else(Region == "Putamen", Channel3_Cells + 295, Channel3_Cells)),
    
    Transduced_Channel3 = if_else(Region == "Caudate", Transduced_Channel3 + 84, 
                                  if_else(Region == "Putamen", Transduced_Channel3 + 79, Transduced_Channel3)),
    
    Channel2_and_3_Cells = if_else(Region == "Caudate", Channel2_and_3_Cells + 49, 
                                   if_else(Region == "Putamen", Channel2_and_3_Cells + 120, Channel2_and_3_Cells)),
    
    Transduced_Channel2_and_3 = if_else(Region == "Caudate", Transduced_Channel2_and_3 + 14, 
                                        if_else(Region == "Putamen", Transduced_Channel2_and_3 + 20, Transduced_Channel2_and_3)),
    
    Total_Neun = if_else(Region == "Caudate", Total_Neun + 0, 
                         if_else(Region == "Putamen", Total_Neun + 0, Total_Neun)),
    
    Neun_Transduced = if_else(Region == "Caudate", Neun_Transduced + 0, 
                              if_else(Region == "Putamen", Neun_Transduced + 0, Neun_Transduced))
  )

joined_summary <- joined_summary %>%
  mutate(
    percentage_transduced_all = Cells_Transduced / Total_Cells * 100,
    percentage_transduced_non_msns = Non_MSNS_Transduced / Total_Non_MSNS * 100,
    percentage_transduced_msns = MSNs_Transduced / Total_MSNs * 100,
    percentage_transduced_ch2 = Transduced_Channel2 / Channel2_Cells * 100,
    percentage_transduced_ch3 = Transduced_Channel3 / Channel3_Cells * 100,
    percentage_transduced_ch2_and_3 = Transduced_Channel2_and_3 / Channel2_and_3_Cells * 100,
    percentage_transduced_neun = ifelse(Total_Neun > 0, Neun_Transduced / Total_Neun * 100, NA),
    .groups = "drop"
  )


# Reshape data for plotting
metric_labels <- c(
  "percentage_transduced_all" = "All Cells",
  "percentage_transduced_non_msns" = "Non-MSN Cells",
  "percentage_transduced_neun" = "Neun Positive Cells",
  "percentage_transduced_msns" = "All MSNs",
  "percentage_transduced_ch2" = "DRD1 Positive MSNs",
  "percentage_transduced_ch3" = "DRD2 Positive MSNs",
  "percentage_transduced_ch2_and_3" = "Double Positive MSNs"
)

# Filter rows with valid values for plotting
joined_summary_long <- joined_summary %>%
  pivot_longer(
    cols = starts_with("percentage"),
    names_to = "Metric",
    values_to = "Percentage"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Percentage = ifelse(is.na(Percentage), 0, Percentage) # Replace NA with 0
  )

individual_summary_long <- individual_summary %>%
  pivot_longer(
    cols = starts_with("percentage"),
    names_to = "Metric",
    values_to = "Percentage"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Percentage = ifelse(is.na(Percentage), 0, Percentage) # Replace NA with 0
  )

individual_summary_long <- individual_summary_long %>%
  filter(!(Animal_ID == "RH16070521" & Metric == "Neun Positive Cells"))

individual_summary_long <- individual_summary_long %>%
  filter(!(Animal_ID == "NA"))

# Define position dodge for bars and points
dodge_position <- position_dodge(width = 0.9)


plot_gp <- ggplot() +
  # Bars: Combined percentages
  geom_bar(
    data = joined_summary_long,
    aes(x = Metric, y = Percentage, fill = Region, group = Region),
    stat = "identity",
    position = dodge_position,
    alpha = 1
  ) +
  # Points: Individual animal percentages
  geom_point(
    data = individual_summary_long,
    aes(x = Metric, y = Percentage, group = Region),
    position = dodge_position,
    size = 2,
    shape = 21,
    color = "black",
    fill = "black",
    show.legend = FALSE
  ) +
  labs(
    title = "Percentage Transduction in GP Animals",
    x = "Cell Type",
    y = "Percent Transduction"
  ) +
  scale_fill_manual(values = c("#6754E2", "#7c7a8c")) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


print(plot_gp)

joined_summary_drop_double <- joined_summary[,-21]
individual_summary_drop_double <- individual_summary[,-22]

# Filter rows with valid values for plotting
joined_summary_long <- joined_summary_drop_double %>%
  pivot_longer(
    cols = starts_with("percentage"),
    names_to = "Metric",
    values_to = "Percentage"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Percentage = ifelse(is.na(Percentage), 0, Percentage) # Replace NA with 0
  )

individual_summary_long <- individual_summary_drop_double %>%
  pivot_longer(
    cols = starts_with("percentage"),
    names_to = "Metric",
    values_to = "Percentage"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Percentage = ifelse(is.na(Percentage), 0, Percentage) # Replace NA with 0
  )

individual_summary_long <- individual_summary_long %>%
  filter(!(Animal_ID == "RH16070521" & Metric == "Neun Positive Cells"))

individual_summary_long <- individual_summary_long %>%
  filter(!(Animal_ID == "NA"))

# Define position dodge for bars and points
dodge_position <- position_dodge(width = 0.9)


plot_no_double <- ggplot() +
  # Bars: Combined percentages
  geom_bar(
    data = joined_summary_long,
    aes(x = Metric, y = Percentage, fill = Region, group = Region),
    stat = "identity",
    position = dodge_position,
    alpha = 1
  ) +
  # Points: Individual animal percentages
  geom_point(
    data = individual_summary_long,
    aes(x = Metric, y = Percentage, group = Region),
    position = dodge_position,
    size = 2,
    shape = 21,
    color = "black",
    fill = "black",
    show.legend = FALSE
  ) +
  labs(
    title = "Percentage Transduction in GP Animals",
    x = "Cell Type",
    y = "Percent Transduction"
  ) +
  scale_fill_manual(values = c("#6754E2", "#7c7a8c")) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(plot_no_double)

individual_summary_long_export <- individual_summary_long[c(1:2,17:18)]
joined_summary_long_export <- joined_summary_long[c(0:1,17:18)]


write.csv(individual_summary_long_export, file = "/Users/christopherfluta/Documents/individual_animals_export.csv")
write.csv(joined_summary_long_export, file = "/Users/christopherfluta/Documents/joined_animals_export.csv")


# Individual plots
individual_summary_long <- individual_summary_long %>%
  filter(!(Animal_ID == "RH16070521"))

Neun_caudate <- individual_summary_long[individual_summary_long$Region == "Caudate", c("Animal_ID", "Region", "Total_Neun", "Neun_Transduced")]
Neun_putamen <- individual_summary_long[individual_summary_long$Region == "Putamen", c("Animal_ID", "Region", "Total_Neun", "Neun_Transduced")]

Neun_caudate <- Neun_caudate[!duplicated(Neun_caudate), ]
Neun_putamen <- Neun_putamen[!duplicated(Neun_putamen), ]

# Calculate total percentage for each region
Neun_caudate$Pct_Transduced <- Neun_caudate$Neun_Transduced / Neun_caudate$Total_Neun * 100
Neun_putamen$Pct_Transduced <- Neun_putamen$Neun_Transduced / Neun_putamen$Total_Neun * 100

Neun_caudate$Total_Pct_Transduced <- sum(Neun_caudate$Neun_Transduced) / sum(Neun_caudate$Total_Neun) * 100
Neun_putamen$Total_Pct_Transduced <- sum(Neun_putamen$Neun_Transduced) / sum(Neun_putamen$Total_Neun) * 100

# Get the unique Total_Pct_Transduced value
total_pct_caudate <- sum(Neun_caudate$Neun_Transduced) / sum(Neun_caudate$Total_Neun) * 100
total_pct_putamen <- sum(Neun_putamen$Neun_Transduced) / sum(Neun_putamen$Total_Neun) * 100

# Create a separate data frame for the bar plot
plot_data_caudate <- data.frame(
  Category = "Caudate",
  Pct_Transduced = total_pct_caudate
)

# Create a separate data frame for the bar plot
plot_data_putamen <- data.frame(
  Category = "Putamen",
  Pct_Transduced = total_pct_putamen
)

# Create the plot
ggplot() +
  # Bar plot for total percentage transduction (ensuring it's only one value)
  geom_bar(data = plot_data_caudate, 
           aes(x = Category, y = total_pct_caudate), 
           stat = "identity", fill = "#6754E2", width = 0.5) +
  
  # Scatter plot for individual animal percentages
  geom_point(data = Neun_caudate, 
             aes(x = "Caudate", y = Pct_Transduced), 
             color = "black", size = 4, alpha = 0.8) +
  
  # Customize labels and theme
  ylim(0, 100) +
  labs(title = "Caudate Neurons Transduced by DB3",
       y = "Percent of Neurons Transduced",
       x = NULL) +
  theme_minimal() +
  theme(text = element_text(family = "IBM Plex Sans", size = 14))

ggsave(
  filename = "Caudate_Neurons_Transduced.pdf", 
  plot = last_plot(),
  device = "pdf"
)

# Create the plot
ggplot() +
  # Bar plot for total percentage transduction (ensuring it's only one value)
  geom_bar(data = plot_data_putamen, 
           aes(x = Category, y = total_pct_putamen), 
           stat = "identity", fill = "#6754E2", width = 0.5) +
  
  # Scatter plot for individual animal percentages
  geom_point(data = Neun_putamen, 
             aes(x = "Putamen", y = Pct_Transduced), 
             color = "black", size = 4, alpha = 0.8) +
  
  # Customize labels and theme
  ylim(0, 100) +
  labs(title = "Putamen Neurons Transduced by DB3",
       y = "Percent of Neurons Transduced",
       x = NULL) +
  theme_minimal() +
  theme(text = element_text(family = "IBM Plex Sans", size = 14))


ggsave(
  filename = "Putamen_Neurons_Transduced.pdf", 
  plot = last_plot(),
  device = "pdf"
)