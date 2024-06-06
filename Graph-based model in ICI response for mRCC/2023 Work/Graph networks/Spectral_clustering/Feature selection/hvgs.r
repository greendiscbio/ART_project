
#-------------------------Import libraries--------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scran")
install.packages("dplyr")

library(scran)
library(dplyr)

#---------------------------Data preparation------------------------------------

# Set working directory
setwd("~/ART_project/2023 Work/Data/Preprocessed_data")

# Load file
file = read.csv('clinic_and_RNA_data_raw_NIVOLUMAB.csv', check.names=FALSE)
print(dim(file))

# Save patients id
patients <- file[,'SUBJID']
print(patients)

# Select only gene variables (122-44015 column)
file <- file[, 123:44015]
cols <- colnames(file)
print(cols[1])
print(cols[43893])
print(dim(file))

# Traspose matrix to insert patients id as name of the column
file <- t(file)
print(dim(file))
colnames(file)<-patients
print(file)

#----------------------Perform High Variable Genes (HVGs)-----------------------

# Not necesary ito normalize since data hhas been already normalized
# Compute log-transformed normalized expression values from a count 
# matrix in a SingleCellExperiment object.
#file <- normalizeCounts(file)
#file

# Model the variance of the log-expression profiles for each gene, 
# decomposing it into technical and biological components based on a 
# fitted mean-variance trend.
dec <- modelGeneVar(file)#, subset.row=1:100)
dec 
dim(dec)

# visualize mean-variance relationship
fit <- metadata(dec)
print(fit)
plot(fit$mean, fit$var, 
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


# Define a set of highly variable genes, based on variance modelling statistics 
# from modelGeneVar or related functions.
top_hvgs <- getTopHVGs(dec, n = 15000)#, prop = 0.4)
print(length(top_hvgs))
print(top_hvgs)

#---------------------------Work with HVGs--------------------------------------
# Create dataframe with top High Variable Genes
file <- t(file)
df_file <- data.frame(file, check.names=FALSE)
print(df_file)
data <- df_file[top_hvgs]
print(dim(data))
data <-t(data)

# visualize mean-variance relationship
dec2 <- modelGeneVar(data)
print(dec2)
fit <- metadata(dec2)
print(fit)
plot(fit$mean, fit$var, 
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


# Compare with Biogrid gene list
data<-t(data)
print(colnames(data))
not_included = c()
all = c()
inclu = c()
setwd("~/ART_project/2023 Work/Graph networks/Biogrid/Data")
rna_genes <- read.table("genes_found.txt", sep = '\n')
rna_genes[["V1"]]
rna_genes <- rna_genes[["V1"]]
for (g in colnames(data)){
  print(g)
  all <- append(all,g)
  included = is.element(g, rna_genes)
  print(included)
  if (included == TRUE){
    inclu <- append(inclu, g)
  }
  else{
    not_included <- append(not_included,g)
  }
}
setwd("~/ART_project/2023 Work/Graph networks/Biogrid/Data/HVGS")
write(not_included, "not_included_genes_hvgs.txt")
write(all, "all_genes_hvgs.txt")
write(inclu, "included_genes_hvgs.txt")
