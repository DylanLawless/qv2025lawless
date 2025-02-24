library(ggplot2)
library(dplyr)

# The following code runs a function on each filtering method and saves the counts of criteria applied.
# We then compare the check of there is complete overlap as validation. 
# Here is the method in more readable terma for just one criteria set:
# source("./guru_signlecase_vcurrent.R")
# source("./ACMG_filters/qv_acmg_filters.R")
# pvs1_table <- table(df$ACMG_PVS1)
# ps1_table <- table(df$ACMG_PS1)
# ps3_table <- table(df$ACMG_PS3)
# ps5_table <- table(df$ACMG_PS5)
# pm2_table <- table(df$ACMG_PM2)
# pm3_table <- table(df$ACMG_PM3)
# 
# psv1_df <- as.data.frame(pvs1_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PVS1")
# ps1_df <- as.data.frame(ps1_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PS1")
# ps3_df <- as.data.frame(ps3_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PS3")
# ps5_df <- as.data.frame(ps5_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PS5")
# pm2_df <- as.data.frame(pm2_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PM2")
# pm3_df <- as.data.frame(pm3_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PM3")
# 
# validation_qv <- bind_rows(psv1_df, ps1_df, ps3_df, ps5_df, pm2_df, pm3_df)
# validation_qv$method <- "QV_yaml"
# print(validation_qv)



# Function to apply ACMG filtering and generate validation data
run_validation <- function(filter_script, method_name) {
  source("./guru_signlecase_vcurrent.R")
  source(filter_script)
  
  # List of ACMG criteria to validate
  acmg_criteria <- c("ACMG_PVS1", "ACMG_PS1", "ACMG_PS3", "ACMG_PS5", "ACMG_PM2", "ACMG_PM3")
  
  # Generate a dataframe for validation counts
  validation_df <- bind_rows(lapply(acmg_criteria, function(criterion) {
    as.data.frame(table(df[[criterion]]), responseName = "Count") %>%
      mutate(Criterion = criterion)
  }))
  
  validation_df$method <- method_name
  return(validation_df)
}

# Run validation for both methods
validation_qv <- run_validation("./ACMG_filters/qv_acmg_filters.R", "QV_yaml")
validation_manual <- run_validation("./ACMG_filters/acmg_filters.R", "QV_manual")

# Combine results
validation_combined <- rbind(validation_qv, validation_manual)

# Print validation results
print(validation_combined)

# Define manual colors
method_colors <- c("QV_yaml" = "#399fe5",  # Blue
                   "QV_manual" = "#ff674b")  # Orange

# Define the desired order of facets
facet_order <- c("ACMG_PVS1", "ACMG_PS1", "ACMG_PS3", "ACMG_PS5", "ACMG_PM2", "ACMG_PM3")
validation_combined$Criterion <- factor(validation_combined$Criterion, levels = facet_order)

p_validation <- ggplot(validation_combined, aes(x = as.factor(Var1), y = Count, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_y_continuous(expand = expansion(mult = 0.8)) + 
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = method_colors) +
  facet_wrap(~ Criterion, scales = "free") +
  labs(x = "ACMG score applied during criteria check", y = "Count of variants lablled", fill = "Method",
       title = "Comparison of ACMG criteria counts on disease cohort",
       subtitle = "QV_acmg_v1.yaml vs Manually encoded") +
  theme_minimal() 

p_validation
print(validation_combined)
print(dim(df))
ggsave(paste(images_directory ,file_suffix, "validation_of_yaml_vs_manual.pdf", sep = "") ,plot = p_validation, width = 8, height = 4)





