# ACMGuru ----
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scico) # devtools::install_github("thomasp85/scico") # scico_palette_show()
library(grid)
library(forcats) # new facet labels
library(ggrepel)
library(ggpubr) # For ggarrange
library(cowplot) # For get_legend
library(patchwork) # plots in panels
library(knitr)

file_suffix <- "Guru_singlecase_"
input_directory <- "../data/"
output_directory <- "../output/"
images_directory <- "../images/"
gnomad_freq_thresh_local <- ""
large_input_directory <- "~/web/spss_exome_vsat/data/ACMGuru_singlecase/"

# acmg ----
# For reference
df_acmg <- fread("../ref/acmg_criteria_table.txt", sep = "\t", header = TRUE, fill=TRUE)
df_acmg_caveat <- fread("../ref/acmg_criteria_table_caveats.txt", sep = "\t", header = TRUE)

# PanelAppRex ----

# Rds format
path_data <- "~/web/PanelAppRex/data/"
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")
rex_core <- readRDS(file= path_PanelAppData_genes_combined_Rds)
colnames(rex_core)[colnames(rex_core) == 'entity_name'] <- 'SYMBOL'
rm(path_data, path_PanelAppData_genes_combined_Rds)

# iuis ----
iuis <- read.table(
  file = "../ref/10875_2022_1289_MOESM2_ESM_DLcleaned.tsv",
  sep = "\t",
  fill = TRUE,  # To handle rows with fewer columns
  header = TRUE # Change this based on whether the first line is a header
)

colnames(iuis)[colnames(iuis) == 'Gene.symbol'] <- 'SYMBOL'

# varsome ----
# LE = less than equal to, GE = greater than equal to
varsome <- read.table(file = "../ref/varsome_calibrated_insilico_thresholds.tsv", sep="\t", header = TRUE)

# qv ----
# Define the chromosome identifiers
chromosomes <- c(1:22, "X")
# chromosomes <- c(21:22, "X") # TEMP TEST
# chromosomes <- c(1)

# Generate file names using paste0 and the chromosome identifiers
file_list <- paste0(
  # input_directory, "bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad_af1_chr_", 
  large_input_directory, "bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad_af1_chr_", 
  chromosomes, 
  ".vcf.gz"
)

df_pathway_list <- list()
# for (f in 21) {
for (f in 1:length(file_list)) {
	cat("Now analysing", f, "\n")
	source("../stand_alone_vcf_to_table/stand_alone_vcf_to_table.R")

	# qv clean ----
	df$cohort_pheno <- df$sample

	# "setpt" = controls "0" / not "setpt" = cases "1"
	df$cohort_pheno[grep("^setpt", df$sample)] <- "0"
	df$cohort_pheno[!grepl("^setpt", df$sample)] <- "1"

	# frequency for cases and controls
	df_genotype_frequency <- df %>%
		dplyr::select(sample, rownames, genotype) %>% 
		unique() %>% # this is import to count genomic positions once rather than transcripts
		mutate(cohort_pheno = ifelse(grepl("^setpt", sample), "0", "1")) %>%
		group_by(rownames, cohort_pheno) %>%
		dplyr::summarize(genotype_total_frequency = sum(genotype)/n(), .groups = "drop") %>%
		pivot_wider(names_from = cohort_pheno, values_from = genotype_total_frequency, names_prefix = "frequency_in_")  %>%
		mutate(is_frequency_in_0_less = ifelse(frequency_in_0 < frequency_in_1, "Yes", "No"))

	df <- df |> filter(genotype > 0) # Keep only variants
	df <- merge(df, df_genotype_frequency, all.x=TRUE)
	rm(df_genotype_frequency)
	
	df <- df |> filter(IMPACT %in% c("HIGH", "MODERATE"))
	
	df <- df |> dplyr::select(-"ClinVar.x",
									  - "ClinVar_CLNSIG.x",
									  - "ClinVar_CLNREVSTAT.x",
									  - "ClinVar_CLNDN.x") # annotation duplicates
	
	df <- df |> distinct()
	df <- df |> filter(cohort_pheno == 1)
	df <- df |> filter(AC < 10)
	
	df_pathway_list[[f]] <- df
}

cat("\nFinished vcf_to_table")

# dfx <- df_pathway_list[[23]] # check chr X

df_pathway <- do.call(rbind, df_pathway_list)
df <- df_pathway
df <- df |> filter(!is.na(SYMBOL)) # clean out unassigned
hold <- df
df <- hold

# rm(list=setdiff(ls(), c("images_directory", "output_directory", "file_suffix", "images_direcory", "rex_core", "gnomad_freq_thresh_local", "df",  "df_acmg", "df_acmg_caveat", "hold", "iuis", "varsome")))
gc()


# rex merge ----
rex_core <- rex_core |> filter(id == "398")
df_rex <- merge(df, rex_core, by="SYMBOL", all.x=TRUE) |> dplyr::select(SYMBOL, everything())

# iuis merge ----
df <- merge(df, iuis, by="SYMBOL", all.x=TRUE) |> dplyr::select(SYMBOL, Inheritance, everything())

# summary ----
# library(Hmisc)
df$gnomAD_AF <- as.numeric(df$gnomAD_AF)
df$AC <- as.numeric(df$AC)
df$AF.x <- as.numeric(df$AF.x)
temp <- df |> ungroup() |> dplyr::select(genotype, Inheritance, IMPACT, Consequence, AF.x, AC, gnomAD_AF, HGVSc) |> unique()

temp |> 
  group_by(genotype) |>
  summarise(n())

temp |> 
  group_by(Inheritance) |>
  summarise(n())

temp |> 
  group_by(IMPACT) |>
  summarise(n())

temp |> 
  group_by(Consequence) |>
  summarise(n())

temp |> 
  group_by(AF.x) |>
  summarise(n())

temp |> 
  group_by(AC) |>
  summarise(n())

temp |>
  ungroup() |>
  dplyr::select(HGVSc) |>
  unique() |>
  summarise(n())

# fix missing chromosome info ----
# chr comes from dbnsfp (I believe), but is a sensible header name so we will reuse it
# Remove "chr" prefix from seqname and create a new 'chr' column
df <- df %>%
  mutate(chr = sub("chr", "", seqnames)) 

# comp_het_flag ----
# flag for comp het. WARNING NOT PHASE CHECKED
df <- df %>%
	group_by(sample, SYMBOL) %>%
	mutate(comp_het_flag = ifelse(n() > 1, 1, NA)) 

# same flag for genotype == 2 (homozygous)
df <- df %>%
	mutate(comp_het_flag = ifelse(is.na(comp_het_flag) & genotype == 2, 1, comp_het_flag)) %>%
	ungroup() %>%
	dplyr::select(comp_het_flag, everything())

# Update the flag for chromosome X observations to include a check for genotype == 2
# NB:: this must be screened using sex info
df <- df %>%
  mutate(comp_het_flag = ifelse(chr == "X" & genotype == 2, 1, comp_het_flag)) %>%
  dplyr::select(comp_het_flag, everything())


# print(n = 30, df |> filter(SYMBOL == "WAS") |> dplyr::select(IMPACT, SYMBOL, genotype))

# source acmg filters ----
print("We are adding PS3 now")

source("./ACMG_filters/distribution_variables.R")
