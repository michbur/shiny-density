library(rentrez)
library(pbapply)
library(dplyr)
library(XML)

# required to properly download xml from ncbi
library(httr)
config(http_version = 2)

# functions to process and annotate data after BEST

if(Sys.info()[["nodename"]] == "amyloid")
  sample1 <- read.csv("/home/michal/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")

gene_names_list <- pblapply(unique(sample1[["prot_id"]])[1L:3], function(ith_name) {

  prot_ids <- entrez_search(db = "protein", term = ith_name)
  config(http_version = 2)
  gene_ids <- entrez_link(dbfrom = "protein", id = prot_ids[["ids"]], db = "gene")
  config(http_version = 0)
  genes <- entrez_fetch(db = "gene", id = gene_ids[["links"]][["protein_gene"]], rettype = "xml")

  gene_list <- xmlToList(genes)
  data.frame(prot_id = ith_name,
             gene_name = gene_list[["Entrezgene"]][["Entrezgene_gene"]][["Gene-ref"]][["Gene-ref_locus"]],
             gene_address = paste0("https://www.ncbi.nlm.nih.gov", 
                                   gene_list[["Entrezgene"]][["Entrezgene_track-info"]][["Gene-track"]][["Gene-track_geneid"]]),
             stringsAsFactors = FALSE)
})
