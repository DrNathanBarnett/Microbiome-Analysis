library(VariantAnnotation)
library(dplyr)
library(data.table)
library(tidyr)
library(caret)
library(keras)
# Read the VCF file
vcf_data <- fread("large_artificial2.vcf", skip = "#CHROM")

# Parse the INFO field into separate columns
vcf_data <- vcf_data %>%
  separate(INFO, into = c("DP", "AF", "GENE", "CONSEQ"), sep = ";") %>%
  mutate(
    DP = as.numeric(sub("DP=", "", DP)),
    AF = as.numeric(sub("AF=", "", AF)),
    GENE = sub("GENE=", "", GENE),
    CONSEQ = sub("CONSEQ=", "", CONSEQ)
  )

# Check the distribution of genes
table(vcf_data$GENE)

# Save to a CSV file for use in the preprocess function
fwrite(vcf_data, "processed_genomic_data.csv")




# Define gene categories
# Epilepsy Genes
epilepsy_genes <- c("ARX", "CDKL5", "SCN1A", "SCN2A", "SCN8A", "STXBP1", "KCNQ2", "KCNT1",
                    "GABRA1", "GABRB3", "GABRG2", "SLC2A1", "SLC6A1", "SLC12A5", "SLC25A22",
                    "GRIN2A", "GRIN2B", "DEPDC5", "ALG13", "PCDH19", "TSC1", "TSC2", "MTOR",
                    "FOXG1", "MECP2", "UBE3A", "SYNGAP1", "CHD2", "DNM1", "KCNB1", "KCNA2",
                    "KCNQ3", "CACNA1A", "CACNA1H", "CACNA1E", "ATP1A2", "ATP1A3", "PIGA",
                    "PNKP", "POLG", "RELN", "TBC1D24", "WDR45",
                    "CACNA2D2", "KCNMA1", "HNRNPU", "SCN3A", "GABRD", "GRIA3")

# Nutrigenomic Genes
nutrigenomic_genes <- c("FADS1", "MTHFR", "VDR", "APOE", "COMT", "GSTM1",
                        "SLC19A1", "SLC25A13", "CBS")

# Pharmacogenomic Genes
pharmacogenomic_genes <- c("CYP2C9", "UGT2B7", "ABCB1", "CYP3A4", "SCN1A", "CYP2D6",
                           "CYP2C19", "CYP1A2", "SLCO1B1", "DPYD")

# Metabolic Genes
metabolic_genes <- c("SLC2A1", "SLC6A1", "ATP1A2", "GLUD1", "ALDH5A1", "LDHB")

# Mitochondrial Genes
mitochondrial_genes <- c("POLG", "PNKP", "NDUFS4", "NDUFV1", "MT-ND5", "SDHA")

# Inflammation Genes
inflammation_genes <- c("IL1B", "IL6", "TNF", "IL10", "CXCL10", "PTGS2")

# Neurotransmitter Genes
neurotransmitter_genes <- c("GAD1", "GAD2", "SLC1A1", "SLC17A7", "SLC32A1", "HTR1A",
                            "HTR2A", "DRD2", "MAOA", "MAOB", "SLC6A3", "SLC6A4")

# Ion Channel Genes
ion_channel_genes <- c("SCN1A", "SCN2A", "SCN8A", "KCNQ2", "KCNT1", "CACNA1A",
                       "CACNA1C", "KCNH2", "KCNN3", "KCNA1")

# Blood-Brain Barrier Genes
bbb_genes <- c("ABCB1", "ABCC1", "SLC2A1", "CLDN5", "OCLN", "VEGFA", "ITGAM",
               "MMP9", "TJP1", "TJP2")

# Pharmaceutical Metabolism Genes
pharmaceutical_metabolism_genes <- c("CYP2C9", "CYP2D6", "CYP3A4", "UGT1A1",
                                     "GSTT1", "SULT1A1", "NAT2", "FMO3")

##########

gene_categories <- list(
  epilepsy = c("ARX", "CDKL5", "SCN1A", "SCN2A", "SCN8A", "STXBP1", "KCNQ2", "KCNT1","GABRA1", "GABRB3", "GABRG2", "SLC2A1", "SLC6A1", "SLC12A5", "SLC25A22","GRIN2A", "GRIN2B", "DEPDC5", "ALG13", "PCDH19", "TSC1", "TSC2", "MTOR","FOXG1", "MECP2", "UBE3A", "SYNGAP1", "CHD2", "DNM1", "KCNB1", "KCNA2","KCNQ3", "CACNA1A", "CACNA1H", "CACNA1E", "ATP1A2", "ATP1A3", "PIGA", "PNKP", "POLG", "RELN", "TBC1D24", "WDR45", "CACNA2D2", "KCNMA1", "HNRNPU","SCN3A", "GABRD", "GRIA3"),
  nutrigenomic = c("FADS1", "MTHFR", "VDR", "APOE", "COMT", "GSTM1", "SLC19A1", "SLC25A13", "CBS"),
  pharmacogenomic = c("CYP2C9", "UGT2B7", "ABCB1", "CYP3A4", "SCN1A", "CYP2D6","CYP2C19", "CYP1A2", "SLCO1B1", "DPYD"),
  metabolic = c("SLC2A1", "SLC6A1", "ATP1A2", "GLUD1", "ALDH5A1", "LDHB"),
  mitochondrial = c("POLG", "PNKP", "NDUFS4", "NDUFV1", "MT-ND5", "SDHA"),
  inflammation = c("IL1B", "IL6", "TNF", "IL10", "CXCL10", "PTGS2"),
  neurotransmitter = c("GAD1", "GAD2", "SLC1A1", "SLC17A7", "SLC32A1", "HTR1A", "HTR2A","DRD2", "MAOA", "MAOB", "SLC6A3", "SLC6A4"),
  ion_channel = c("SCN1A", "SCN2A", "SCN8A", "KCNQ2", "KCNT1", "CACNA1A", "CACNA1C","KCNH2", "KCNN3", "KCNA1"),
  bbb = c("ABCB1", "ABCC1", "SLC2A1", "CLDN5", "OCLN", "VEGFA", "ITGAM", "MMP9", "TJP1", "TJP2"),
  pharmaceutical_metabolism = c("CYP2C9", "CYP2D6", "CYP3A4", "UGT1A1","GSTT1", "SULT1A1", "NAT2", "FMO3")
)

# Filter genomic data for each gene category
epilepsy_data <- filter_genes(genomic_data, epilepsy_genes, dataset_type = "epilepsy")
nutrigenomic_data <- filter_genes(genomic_data, nutrigenomic_genes, dataset_type = "nutrigenomic")
pharmacogenomic_data <- filter_genes(genomic_data, pharmacogenomic_genes, dataset_type = "pharmacogenomic")

# Combine all filtered datasets
combined_data <- combine_filtered_data(epilepsy_data, nutrigenomic_data, pharmacogenomic_data)




genomic_file <- "Path"
# Preprocess Genomic Data
preprocess_genomic_data <- function(genomic_file) {
  genomic_data <- fread(genomic_file)
  genomic_data <- genomic_data %>%
    mutate(impact_score = case_when(
      CONSEQ %in% c("missense_variant", "stop_gained") ~ 3,
      CONSEQ %in% c("splice_region_variant", "inframe_insertion") ~ 2,
      TRUE ~ 1
    ))
  return(genomic_data)
}

#genomic_data <- preprocess_genomic_data(genomic_file)

# Filter Data for Specific Gene Categories

filter_genes <- function(genomic_data, gene_list, dataset_type) {
  genomic_data %>%
    filter(GENE %in% gene_list) %>%
    mutate(
      impact_score = case_when(
        CONSEQ %in% c("missense_variant", "stop_gained") ~ sample(2:4, n(), replace = TRUE), # Higher impact for critical variants
        CONSEQ %in% c("splice_region_variant", "inframe_insertion") ~ sample(1:3, n(), replace = TRUE), # Medium impact
        TRUE ~ sample(0:2, n(), replace = TRUE) # Lower impact
      ),
      dataset_type = dataset_type # Add dataset type for traceability
    ) %>%
    group_by(SAMPLEID, dataset_type) %>%
    summarize(
      total_impact_score = sum(impact_score),
      variant_count = n(),
      .groups = "drop"
    )
}

filtered_data <- lapply(names(gene_categories), function(category) {
  filter_genes(genomic_data, gene_categories[[category]], category)
})

# Combine all filtered data
combine_ddata <- do.call(rbind, filtered_data)

# Combine Filtered Data
#combine_filtered_data <- function(...) {
#datasets <- list(...) # Capture multiple datasets as a list
#combined_data <- datasets %>%
#purrr::reduce(function(x, y) {
# Ensure SAMPLEID exists in both datasets
#if (!"SAMPLEID" %in% colnames(x)) x <- x %>% mutate(SAMPLEID = paste0("Sample", row_number()))
#if (!"SAMPLEID" %in% colnames(y)) y <- y %>% mutate(SAMPLEID = paste0("Sample", row_number()))

# Merge datasets
#left_join(x, y, by = "SAMPLEID", suffix = c("_x", "_y"))
#}) %>%
#mutate(across(where(is.numeric), ~ coalesce(.x, 0))) # Fill NAs with 0 for numeric columns
#return(combined_data)
#}
# Feature Engineering

colnames(combined_data)

feature_engineering <- function(combined_data, clinical_data) {
  features <- combined_data %>%
    left_join(clinical_data, by = "SAMPLEID") %>%
    mutate(
      # Calculate risk scores separately for each category
      risk_score_epilepsy = total_impact_score_x,
      risk_score_nutrigenomic = total_impact_score_y,
      risk_score_pharmacogenomic = total_impact_score,
      
      # Assign risk categories for each category
      risk_alert_epilepsy = case_when(
        risk_score_epilepsy > threshold ~ "High Risk",
        risk_score_epilepsy > threshold / 2 ~ "Moderate Risk",
        risk_score_epilepsy > threshold / 4 ~ "Low Risk",
        TRUE ~ "No Risk"
      ),
      risk_alert_nutrigenomic = case_when(
        risk_score_nutrigenomic > threshold ~ "High Risk",
        risk_score_nutrigenomic > threshold / 2 ~ "Moderate Risk",
        risk_score_nutrigenomic > threshold / 4 ~ "Low Risk",
        TRUE ~ "No Risk"
      ),
      risk_alert_pharmacogenomic = case_when(
        risk_score_pharmacogenomic > threshold ~ "High Risk",
        risk_score_pharmacogenomic > threshold / 2 ~ "Moderate Risk",
        risk_score_pharmacogenomic > threshold / 4 ~ "Low Risk",
        TRUE ~ "No Risk"
      )
    )
  return(features)
}

#features <- feature_engineering(combined_data, clinical_data)

# Save the data into an RDS file for use in the RMarkdown report
saveRDS(list(
  filtered_data = combine_ddata,
  quality_metrics = quality_metrics,
  impact_analysis = impact_analysis,
  feature_data = features
), "report_data.rds")"