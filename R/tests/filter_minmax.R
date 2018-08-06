app <- ShinyDriver$new("../")
app$snapshotInit("filter_minmax")

app$setInputs(method = "Minimum/maximum values",
              var = "mal_WT",
              type = "2",
              x = 8,
              filter = "click")

app$snapshot()
