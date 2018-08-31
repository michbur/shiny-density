library(rentrez)
library(pbapply)
library(dplyr)
library(XML)

# required to properly download xml from ncbi
httr::set_config(httr::config(http_version = 0))

# functions to process and annotate data after BEST

if(Sys.info()[["nodename"]] %in% c("amyloid", "lori"))
  sample1 <- read.csv("/home/michal/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")

gene_names_list <- pblapply(unique(sample1[["prot_id"]]), function(ith_name) try({

  prot_ids <- entrez_search(db = "protein", term = ith_name)
  #config(http_version = 2)
  gene_ids <- entrez_link(dbfrom = "protein", id = prot_ids[["ids"]], db = "gene")
  #config(http_version = 0)
  genes <- entrez_fetch(db = "gene", id = gene_ids[["links"]][["protein_gene"]], rettype = "xml")

  gene_list <- xmlToList(genes)
  data.frame(prot_id = ith_name,
             gene_name = gene_list[["Entrezgene"]][["Entrezgene_gene"]][["Gene-ref"]][["Gene-ref_locus"]],
             gene_address = paste0("https://www.ncbi.nlm.nih.gov/", 
                                   gene_list[["Entrezgene"]][["Entrezgene_track-info"]][["Gene-track"]][["Gene-track_geneid"]]),
             stringsAsFactors = FALSE)
}, silent = TRUE))


if(Sys.info()[["nodename"]] %in% c("amyloid", "lori"))
  save(gene_names_list, file = "/home/michal/Dropbox/PepArray_results/2018-06-07/gene_names_list.RData")

if(Sys.info()[["nodename"]] %in% c("amyloid", "lori"))
  load("/home/michal/Dropbox/PepArray_results/2018-06-07/gene_names_list.RData")

gene_names_df <- gene_names_list[!sapply(gene_names_list, class) == "try-error"] %>% 
  bind_rows() %>% 
  mutate(gene_name = paste0('<a href="', gene_address, 
                            '" target="_blank" class="btn btn-primary">', gene_name, '</a>')) %>% 
  select(-gene_address)


# final object ----------------------------------------

full_data <- left_join(sample1, gene_names_df)

if(Sys.info()[["nodename"]] %in% c("amyloid", "lori"))
  save(full_data, file = "/home/michal/Dropbox/PepArray_results/2018-06-07/full_data.RData")
