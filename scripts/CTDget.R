# Utility for converting CTD to ksRepo list format
# Author: Adam Brown
# Date: 12/2/15

# Download gene-compound interaction list from http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz
# Remove header from CTD_chem_gene_ixns.csv before running this script

# Read csv
CTD <- read.table(file = 'CTD_chem_gene_ixns.csv', sep = ',', header = T, quote='\"', stringsAsFactors = F)
CTD.human <- subset(CTD, Organism == 'Homo sapiens') # Only human interactions
comps <- unique(CTD.human$ChemicalName) # Get unique compounds

# Produce ksRepo compatible list
CTD.list <- vector(mode = 'list', length = length(comps)) # Initialize
names(CTD.list) <- comps # Initialize

for (i in 1:length(comps)) {
  CTD.list[[i]] <- CTD$GeneID[which(CTD$ChemicalName == comps[i])] # Make output
}

CTD.list <- lapply(CTD.list, unique) # Trim