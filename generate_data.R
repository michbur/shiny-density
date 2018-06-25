library(dplyr)

set.seed(15390)

dat <- data.frame(seq = sapply(1L:1000, function(dummy) paste0(sample(LETTERS, 10, replace = TRUE), collapse = "")),
           ben_mal = rnorm(1L:1000),
           ben_WT = rnorm(1L:1000),
           mal_WT = rnorm(1L:1000))

reshape2::melt(dat, value.name = "InROPE", variable.name = "type") %>% 
  write.csv("./data/sample1.csv", row.names = FALSE)

