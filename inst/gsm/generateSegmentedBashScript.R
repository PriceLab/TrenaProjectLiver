max <- 13744
starts <- seq(1, 13760, 40)
ends <- starts + 39
ends[length(ends)] <- max
bashScript <- file("runByChunks.sh")
lines <- sprintf("Rscript stagedRun.R %d %d", starts, ends)
writeLines(lines, bashScript)
close(bashScript)
