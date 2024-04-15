library(data.table)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data
Bac <- fread("Data/otu_table_bacteria.txt")
Fung <- fread("Data/otu_table_fungi.txt")

# Sample names
Bac_names <- colnames(Bac)[-1]
Fung_names <- colnames(Fung)[-1]

# Adjusting sample names of ITS data to harmonize with 16S data
num <- paste0("_", 1:6)
num2 <- as.character(1:6)
for(i in 1:6){
  Fung_names[grepl(num[i], Fung_names)] <- gsub(num[i],num2[i],Fung_names[grepl(num[i], Fung_names)])
}

# Samples common to both datasets
Common_names <- Bac_names[Bac_names %in% Fung_names]
# Experiment number
Exp <- as.numeric(substr(Common_names, 1, 1))
# Compartment
Compartment <- fcase(nchar(Common_names) == 4, "Soil",
                     nchar(Common_names) == 5, "Root",
                     nchar(Common_names) == 6, "Rhizosphere")
# Genotype
Genotype <- fcase(substr(Common_names,3,3) == "G", "Gifu",
                  substr(Common_names,3,3) == "N", "nfr5",
                  substr(Common_names,3,3) == "R", "ram1",
                  substr(Common_names,3,3) == "S", "symrk",
                  substr(Common_names,3,3) == "C", "ccamk")
# Samples with four-character names are soil samples
Genotype[nchar(Common_names) == 4] <- "Soil"

# Collecting in data table
Meta_data <- data.table(SampleID = Common_names,
                        Experiment = Exp,
                        Compartment = Compartment,
                        Genotype = Genotype)
fwrite(Meta_data, "Data/meta_data.csv")
