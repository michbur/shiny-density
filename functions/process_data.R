library(rentrez)
library(pbapply)
library(dplyr)
library(XML)

# functions to process and annotate data after BEST

if(Sys.info()[["nodename"]] == "amyloid")
  sample1 <- read.csv("/home/michal/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")

gene_names_list <- pblapply(unique(sample1[["prot_id"]]), function(ith_name) {
  prot_ids <- entrez_search(db = "protein", term = ith_name)
  gene_ids <- entrez_link(dbfrom = "protein", id = prot_ids[["ids"]], db = "gene")
  genes <- entrez_fetch(db = "gene", id = gene_ids[["links"]][["protein_gene"]], rettype = "fasta")

  gene_lines <- strsplit(genes, "\n")[[1]]
  
  data.frame(prot_id = ith_name,
             gene_name = strsplit(gene_lines[2], " ")[[1]][2],
             gene_address = paste0("https://www.ncbi.nlm.nih.gov/", 
                                   gsub("ID: ", "", gene_lines[grep(pattern = "ID: ", x = gene_lines)])),
             stringsAsFactors = FALSE)
}) 
