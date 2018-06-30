library(dplyr)

set.seed(15390)

dat <- data.frame(seq = sapply(1L:1000, function(dummy) paste0(sample(LETTERS, 10, replace = TRUE), collapse = "")),
                  ben_mal = rnorm(1L:1000),
                  ben_WT = rnorm(1L:1000),
                  mal_WT = rnorm(1L:1000),
                  prot_id = sample(paste0("prot", 1L:200), size = 1000, replace = TRUE))

reshape2::melt(dat, value.name = "InROPE", variable.name = "type") %>% 
  write.csv("./data/sample1.csv", row.names = FALSE)

