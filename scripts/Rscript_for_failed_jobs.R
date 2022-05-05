# Command to extract JobIDs and status for a given jobID
# sacct -n -X -j 44417975 -o JobID%-13,state | sort | uniq -c > 4417975_failed_jobs.txt

# args[1] is the txt file to laod
args = commandArgs(trailingOnly=TRUE)

jobs <- read.csv(args[1], header = F, sep = " ")

jobs <- jobs[,c("V8", "V10")]

jobs$V11 <- substr(jobs$V8, start = 10, stop = 15)

failed_ids <- jobs$V11[jobs$V10 != "COMPLETED"]
failed_ids <- failed_ids[failed_ids != ""]
failed_ids <- sort(failed_ids)
write.table(failed_ids[3:length(failed_ids)], file = paste("failed_ids.csv", args[1], sep = "_"), quote = F, col.names = F, sep = ",")

print(length(failed_ids))
print(noquote(paste(failed_ids, collapse = ",")))
