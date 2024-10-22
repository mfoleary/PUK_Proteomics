# User Instructions: 
# Replace the blank file paths below with the appropriate paths on your system.

# Set working directory
setwd("") 

# Define the file paths (update with the correct file locations)
data_dir <- ""         # Folder containing data files
data_file <- ""        # CSV file for gene expression matrix
metadata_file <- ""    # CSV file for metadata
plots_dir <- ""        # Folder to save plots
results_dir <- ""      # Folder to save results
cor_file <- ""         # CSV file to save module-trait correlations
pvals_file <- ""       # CSV file to save p-values of module-trait correlations
output_file <- ""      # Excel file to save filtered data
Metacor_file <- ""     # CSV file to save metadata correlations
Metapvals_file <- ""   # CSV file to save metadata correlation p-values

# Create a user library directory (change the path as needed)
user_library <- ""

# Install necessary packages if not already installed
required_packages <- c("BiocManager", "GO.db", "DESeq2", "WGCNA", "ggforce", "ggplot2", 
                       "magrittr", "readr", "dplyr", "stringr", "tibble", "openxlsx", "limma")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# For Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GO.db", "DESeq2", "limma"))


library(BiocManager)
library(GO.db) 
library(DESeq2)
library(WGCNA)
library(ggforce)
# We'll be making some plots
library(ggplot2)
# We will need this so we can use the pipe: %>%
library(magrittr)
library(readr)

# Define the file path to the data directory
# Replace with the path of the folder the files will be in
data_dir <-""

# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
# Replace with the path to your dataset file

data_file <- ""
metadata_file <- ""

# Define the file path to the plots directory
plots_dir <- ""

# Define path to results
results_dir <- ""

# Read in metadata CSV file
metadata <- readr::read_csv(metadata_file)
# Read in data CSV file
df <- readr::read_csv(data_file) %>%
  # Here we are going to store the gene IDs as row names so that we can have a numeric matrix to perform calculations on later
  tibble::column_to_rownames("#Uniprot")

#Needed for next step
library(dplyr)  

# Make the data in the order of the metadata. You cannot have a space in the metadata file column heading. 
df <- df %>%
  dplyr::select(metadata$Participant)

# Check if this is in the same order
all.equal(colnames(df), metadata$Participant)

# Round the values in the data frame. This is why is it important to import raw proteomics numbers and transform later
df <- round(df)

metadata <- metadata %>%
  dplyr::mutate(
    Timepoint = dplyr::case_when(
      
      stringr::str_detect(Participant, "B") ~ "Baseline",
      stringr::str_detect(Participant, "T") ~ "Trained",
    ),
    
    # Convert to factor
    Timepoint = as.factor(Timepoint)
  )

levels(metadata$Timepoint)

#Load DESeq2
library(DESeq2)

# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = df, # Our prepped data frame with counts
  colData = metadata, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)

# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)

# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data

#This step might take a while, don't panic
sft <- pickSoftThreshold(normalized_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()


bwnet <- blockwiseModules(normalized_counts,
                          maxBlockSize = 5000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = 10, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          # so we should set a seed
)

readr::write_rds(bwnet,
                 file = file.path("")
)

module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

all.equal(metadata$Participant, rownames(module_eigengenes))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(BiocManager)
library(limma)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

head(stats_df)

gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

# Print out a preview
head(module_eigengenes)

head(colors)

#Timepoint Bin

# Convert 'Timepoint' column to character
metadata <- metadata %>%
  mutate(Timepoint = as.character(Timepoint))

metadata <- metadata %>%
  mutate(Timepoint_bin = ifelse(grepl('Trained', Timepoint, ignore.case = TRUE), 1, 0))

head(metadata)

metadata <- metadata %>%
  mutate(Timepoint_bin = ifelse(grepl('Trained', Timepoint), 1, 0))

# Create 'Timepoint_bin' column based on 'Timepoint' content
metadata <- metadata %>%
  mutate(Timepoint_bin = case_when(
    grepl('Trained', Timepoint, ignore.case = TRUE) ~ 1,
    TRUE ~ 0  # Assign 0 for all other cases
  ))

nsamples <- nrow(metadata)
nGenes <- nrow(df)

print(nsamples)
print(nGenes)

traits <- metadata 

library(tibble)

# Assuming "Participant" is a column in your traits tibble
traits <- traits %>%
  column_to_rownames(var = "Participant")

module.trait.cor <- cor(module_eigengenes, traits, use = 'p')
module.trait.cor.pvals<- corPvalueStudent(module.trait.cor, nsamples)

# Define the file paths
cor_file <- ""
pvals_file <- ""

# Write module.trait.cor to CSV
write.csv(module.trait.cor, file = cor_file, row.names = FALSE)

# Write module.trait.cor.pvals to CSV
write.csv(module.trait.cor.pvals, file = pvals_file, row.names = FALSE)

gene_module_key %>%
  dplyr::filter(module == "ME12")
gene_module_key %>%
  dplyr::filter(module == "ME11")
gene_module_key %>%
  dplyr::filter(module == "ME7")
gene_module_key %>%
  dplyr::filter(module == "ME6")
gene_module_key %>%
  dplyr::filter(module == "ME10")
gene_module_key %>%
  dplyr::filter(module == "ME9")
gene_module_key %>%
  dplyr::filter(module == "ME1")
gene_module_key %>%
  dplyr::filter(module == "ME8")
gene_module_key %>%
  dplyr::filter(module == "ME3")
gene_module_key %>%
  dplyr::filter(module == "ME5")
gene_module_key %>%
  dplyr::filter(module == "ME2")
gene_module_key %>%
  dplyr::filter(module == "ME4")
gene_module_key %>%
  dplyr::filter(module == "ME0")

filtered_data <- gene_module_key %>%
  dplyr::filter(module == "ME0")

print(filtered_data)
warning


library(dplyr)
library(openxlsx)

# Define the file path
output_file <- ""

# Filter data for each module
module_ME12 <- gene_module_key %>%
  filter(module == "ME12")

module_ME11 <- gene_module_key %>%
  filter(module == "ME11")

module_ME7 <- gene_module_key %>%
  filter(module == "ME7")

module_ME6 <- gene_module_key %>%
  filter(module == "ME6")

module_ME10 <- gene_module_key %>%
  filter(module == "ME10")

module_ME9 <- gene_module_key %>%
  filter(module == "ME9")

module_ME1 <- gene_module_key %>%
  filter(module == "ME1")

module_ME8 <- gene_module_key %>%
  filter(module == "ME8")

module_ME3 <- gene_module_key %>%
  filter(module == "ME3")

module_ME5 <- gene_module_key %>%
  filter(module == "ME5")

module_ME2 <- gene_module_key %>%
  filter(module == "ME2")

module_ME4 <- gene_module_key %>%
  filter(module == "ME4")

module_ME0 <- gene_module_key %>%
  filter(module == "ME0")

# Write to Excel file
write.xlsx(
  list(ME12 = module_ME12, ME11 = module_ME11, ME7 = module_ME7, ME6 = module_ME6, ME10 = module_ME10,
       ME9 = module_ME9, ME1 = module_ME1, ME8 = module_ME8, ME3 = module_ME3, ME5 = module_ME5,
       ME2 = module_ME2, ME4 = module_ME4, ME0 = module_ME0),
  file = output_file,
  sheetName = c("ME12", "ME11", "ME7", "ME6", "ME10",
                "ME9", "ME1", "ME8", "ME3", "ME5",
                "ME2", "ME4", "ME0"),
  row.names = FALSE
)


#Compute Pearson correlations between metadata traits (not WGCNA modules)
trait_cor <- cor(traits, use = 'p')

#Calculate p-values for the correlations
library(WGCNA)
trait_cor_pvals <- corPvalueStudent(trait_cor, 28)

# Define the file paths
Metacor_file <- ""
Metapvals_file <- ""

# Write module.trait.cor to CSV
write.csv(trait_cor, file = Metacor_file, row.names = FALSE)

# Write module.trait.cor.pvals to CSV
write.csv(trait_cor_pvals, file = Metapvals_file, row.names = FALSE)

