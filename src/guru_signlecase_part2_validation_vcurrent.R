library(ggplot2)
library(dplyr)

# validate two methods ----
## method 1 ----
source("./guru_signlecase_vcurrent.R")
source("./ACMG_filters/qv_acmg_filters.R")
ps1_table <- table(df$ACMG_PS1)
ps3_table <- table(df$ACMG_PS3)
ps5_table <- table(df$ACMG_PS5)
pm2_table <- table(df$ACMG_PM2)

ps1_df <- as.data.frame(ps1_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PS1")
ps3_df <- as.data.frame(ps3_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PS3")
ps5_df <- as.data.frame(ps5_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PS5")
pm2_df <- as.data.frame(pm2_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PM2")

validation_qv <- bind_rows(ps1_df, ps3_df, ps5_df, pm2_df)
validation_qv <- bind_rows(ps1_df, ps3_df, ps5_df)
validation_qv$method <- "QV_yaml"
print(validation_qv)

## method 2 ----
source("./guru_signlecase_vcurrent.R")
source("./ACMG_filters/acmg_filters.R")
ps1_table <- table(df$ACMG_PS1)
ps3_table <- table(df$ACMG_PS3)
ps5_table <- table(df$ACMG_PS5)
pm2_table <- table(df$ACMG_PM2)

ps1_df <- as.data.frame(ps1_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PS1")
ps3_df <- as.data.frame(ps3_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PS3")
ps5_df <- as.data.frame(ps5_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PS5")
pm2_df <- as.data.frame(pm2_table, responseName = "Count") %>% mutate(Criterion = "ACMG_PM2")

validation_manual <- bind_rows(ps1_df, ps3_df, ps5_df, pm2_df)
validation_manual$method <- "QV_manual"
print(validation_manual)

# Compare -----
validation_combined <- rbind(validation_qv, validation_manual)

# Define manual colors
method_colors <- c("QV_yaml" = "#399fe5",  # Blue
                   "QV_manual" = "#ff674b")  # Orange

p_validation <- ggplot(validation_combined, aes(x = as.factor(Var1), y = Count, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = method_colors) +
  facet_wrap(~ Criterion, scales = "free") +
  labs(x = "ACMG score applied during criteria check", y = "Count of variants lablled", fill = "Method",
       title = "Comparison of ACMG criteria counts on disease cohort",
       subtitle = "QV_acmg_v1.yaml vs Manually encoded") +
  theme_minimal() 

p_validation
print(validation_combined)
print(dim(df))
ggsave(paste(images_directory ,file_suffix, "validation_of_yaml_vs_manual.pdf", sep = "") ,plot = p_validation, width = 8, height = 3)

# > print(validation_combined)
# Var1 Count Criterion    method
# 1     0   376  ACMG_PS1   QV_yaml
# 2     4    47  ACMG_PS1   QV_yaml
# 3     0   308  ACMG_PS3   QV_yaml
# 4     4   115  ACMG_PS3   QV_yaml
# 5     0   419  ACMG_PS5   QV_yaml
# 6     4     4  ACMG_PS5   QV_yaml
# 7     0   376  ACMG_PS1 QV_manual
# 8     4    47  ACMG_PS1 QV_manual
# 9     0   308  ACMG_PS3 QV_manual
# 10    4   115  ACMG_PS3 QV_manual
# 11    0   419  ACMG_PS5 QV_manual
# 12    4     4  ACMG_PS5 QV_manual
# > ggsave(paste(images_directory ,file_suffix, "validation_of_yaml_vs_manual.pdf", sep = "") ,plot = p_validation, width = 8, height = 3)
# > print(dim(df))
# [1] 423 377
