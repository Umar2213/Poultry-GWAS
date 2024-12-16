# Poultry-GWAS
Genome-Wide Association Study for Poultry A comprehensive analysis pipeline including SNP interactions and predictive modeling
R-based GWAS Pipeline:
1. Loading and Preprocessing the Data
library(data.table)  # For fast data manipulation
library(dplyr)       # For data manipulation
library(gdsfmt)      # For Genomic Data Structure (GDS) file format
library(SNPRelate)   # For GWAS analysis
library(qqman)       # For plotting GWAS results
phenotype_data <- fread("poultry_phenotypes.csv")
genotype_data <- snpgdsOpen("poultry_genotypes.gds")
phenotype_data <- phenotype_data %>%
  filter(ID %in% snpgdsSampleID(genotype_data))  # Match IDs between phenotypic and genotypic data
2. GWAS Analysis - Basic SNP Association
result <- snpgdsGWAS(
  genofile = genotype_data, 
  phenofile = phenotype_data$Trait, 
  method = "linear", 
  covar = NULL,  # You can include covariates (e.g., age, sex) if applicable
  threads = 4    # Parallel processing
)
significant_SNPs <- result$results %>%
  filter(P.value < 5e-8)  # GWAS threshold for significance
write.csv(significant_SNPs, "significant_SNPs_poultry.csv")
3. Custom Metric - SNP Interaction Analysis 

interaction_analysis <- function(genotype_data, phenotype_data) {
  snp_pairs <- combn(snpgdsSNPRate(genotype_data), 2, simplify = TRUE)
  interaction_pvalues <- c()
  for (pair in 1:ncol(snp_pairs)) {
    snp1 <- snp_pairs[1, pair]
    snp2 <- snp_pairs[2, pair]
    interaction_model <- lm(phenotype_data$Trait ~ genotype_data[, snp1] * genotype_data[, snp2])
    interaction_pvalues[pair] <- summary(interaction_model)$coefficients[3, 4]  # P-value for interaction term
  }
  return(data.frame(snp1 = snp_pairs[1, ], snp2 = snp_pairs[2, ], pvalue = interaction_pvalues))
}
interaction_results <- interaction_analysis(genotype_data, phenotype_data)
write.csv(interaction_results, "SNP_interactions_poultry.csv")






4. Unique Plotting Method for Poultry GWAS - Trait Specific
ggplot(result$results, aes(x = Position, y = -log10(P.value))) +
  geom_point(aes(color = factor(Chromosome)), size = 1) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "GWAS Results for Poultry Trait", 
       x = "Genomic Position", 
       y = "-log10(P-value)",
       color = "Chromosome") +
  scale_color_manual(values = c("blue", "green", "orange", "purple", "yellow"))  # Custom colors for chromosomes
5. SNP Prioritization - Using Gene Annotation
gene_annotations <- fread("gene_annotations.csv")
find_nearest_gene <- function(snp_position, gene_annotations) {
  distances <- abs(gene_annotations$start_position - snp_position)
  closest_gene <- gene_annotations[which.min(distances), ]
  return(closest_gene$gene_name)
}
significant_SNPs$nearest_gene <- sapply(significant_SNPs$Position, find_nearest_gene, gene_annotations)
write.csv(significant_SNPs, "prioritized_SNPs_poultry.csv")
6. Cross-Validation for Model Performance 
trainIndex <- createDataPartition(phenotype_data$Trait, p = 0.8, list = FALSE)
train_data <- phenotype_data[trainIndex, ]
test_data <- phenotype_data[-trainIndex, ]
predictive_model <- lm(Trait ~ ., data = train_data[, c("Trait", significant_SNPs$ID)])
predictions <- predict(predictive_model, test_data)
mse <- mean((predictions - test_data$Trait)^2)
rsq <- summary(predictive_model)$r.squared
cat("Mean Squared Error:", mse, "\n")
cat("R-squared:", rsq, "\n")
![image](https://github.com/user-attachments/assets/1b496feb-1e1a-451f-bd55-592fe4c81fed)
