
# This is a few lines to make the index_array for either 40x or 1000x intervals and whichever samples I'm using (just 5 or all)

# intervals <- list.files("~/projects/cichlids_2/genome/intervals_ASM185804v2_1000x/")
samples <- list.files("~/ncbi/public/sra/")
# test <- lapply(samples, function(x) data.frame(sample = x, interval = intervals))
# test <- lapply(samples[1:5], function(x) data.frame(sample = x, interval = intervals))
# test2 <- Reduce(rbind, test)

write.table(samples, "~/projects/cichlids_2/scripts/index_samples.csv", quote = FALSE, col.names = FALSE, sep = ",")

