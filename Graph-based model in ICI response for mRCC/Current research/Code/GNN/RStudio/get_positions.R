if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BioMart
# BiocManager::install("biomaRt")

library(biomaRt)
library("xlsx")

# Create a query in BioMart with the genes of interest
hsapiens37_mart <- useEnsembl('ensembl', GRCh = 37, dataset = 'hsapiens_gene_ensembl') # important the ref genome! GRCh37

# Output attributes
#out_attributes <- c('external_gene_name', 'external_transcript_name','ensembl_gene_id', 'ensembl_transcript_id', 'description', 'chromosome_name', 'start_position', 'end_position', 'transcript_start', 'transcript_end')
out_attributes <- c('external_gene_name','ensembl_gene_id', 'description', 'chromosome_name', 'start_position', 'end_position' )


################################################################################
# Obtain BED files for AD gene set
################################################################################

infile <- 'mRCc_genes_disgenet.txt'

# Filters - inputs for a biomaRt query
gene_filter <- c('external_gene_name')
gene_names <- scan(infile, what = '', sep = '\n')

filters = listAttributes(hsapiens37_mart)
filters[1:25,]

# getBM - main biomaRt query function
out_df <- getBM(attributes = out_attributes,
                filters = gene_filter,
                values = gene_names,
                mart = hsapiens37_mart)
out_df
write.xlsx(out_df,
            col.names = TRUE,
            row.names = FALSE,
            file = 'mapping_results.xlsx')

