set.seed(15390)

dat <- data.frame(seq = sapply(1L:1000, function(dummy) paste0(sample(LETTERS, 10, replace = TRUE), collapse = "")),
           ben_mal = rnorm(1L:1000),
           ben_WT = rnorm(1L:1000),
           mal_WT = rnorm(1L:1000))

write.csv(dat, "./data/sample1.csv", row.names = FALSE)

