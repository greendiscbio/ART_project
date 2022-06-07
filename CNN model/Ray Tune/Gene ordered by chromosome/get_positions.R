#! /usr/bin/Rscript

#setwd('/home/laura/Desktop')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
# BioMart
library(biomaRt)

# Create a query in BioMart with the genes of interest
hsapiens37_mart <- useEnsembl('ensembl', GRCh = 37, dataset = 'hsapiens_gene_ensembl') # important the ref genome! GRCh37

# Output attributes
out_attributes <- c('external_gene_name','ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position')


################################################################################
# Obtain BED files for AD gene set
################################################################################

infile <- 'symbol.txt'

# Filters - inputs for a biomaRt query
gene_filter <- c('external_gene_name')
gene_names <- scan(infile, what = '', sep = '\n')

# getBM - main biomaRt query function
out_df <- getBM(attributes = out_attributes,
                filters = gene_filter,
                values = gene_names,
                mart = hsapiens37_mart)
out_df
write.table(out_df,
            sep = '\t',
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            file = 'positions2.tsv')




