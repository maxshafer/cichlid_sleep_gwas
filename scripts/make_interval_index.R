
# This is a few lines to make the index_array for either 40x or 1000x intervals and whichever samples I'm using (just 5 or all)

intervals <- list.files("~/projects/cichlids_2/genome/intervals_ASM185804v2_chr/")
samples <- list.files("~/ncbi/public/sra/")
test <- lapply(samples, function(x) data.frame(sample = x, interval = intervals))
test <- lapply(samples, function(x) data.frame(sample = x, interval = intervals))
test2 <- Reduce(rbind, test)

write.table(test2, "~/projects/cichlids_2/scripts/index_array_chr.csv", quote = FALSE, col.names = FALSE, sep = ",")

