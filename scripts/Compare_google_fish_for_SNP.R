library(gsheet)
library(seqinr)
library(vcd)

# Set the working directory (where I download the fasta files)
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/")

# Can use this, but maybe doesn't work as well as raw data
resolved_names <- read.csv("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Fish_sleep/resolved_names_local.csv")

# Load the google sheet
url <- 'https://docs.google.com/spreadsheets/d/18aNqHT73hX06cGRlf6oj7Y4TVKf6jd_Q5ojNIxm2rys/edit#gid=0'
sleepy_fish <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
sleepy_fish$Diel_Pattern <- tolower(sleepy_fish$Diel_Pattern)
sleepy_fish <- sleepy_fish[sleepy_fish$Diel_Pattern != "",]


# Load the fast files (list of fastas)

fasta.file.names <- list.files()[grep("alignment", list.files())]

fasta.file <- read.fasta("LOC100696920_G659R_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("LOC109204088_T51F_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("susd1_G895R_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("cep68_G360A_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("LOC102077138_N406H_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("sim2_A610T_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("LOC100705676_V339I_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("pkhd1l1_K3099R_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("col4a1_D1252G_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("LOC100696545_S677L_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("LOC100695019_I32F_alignment", seqtype = "AA", set.attributes = T, whole.header = T)
fasta.file <- read.fasta("LOC100701194_A331V_alignment", seqtype = "AA", set.attributes = T, whole.header = T)

# fasta.file[[1]][1] # this is the amino acid
# attr(fasta.file[[1]], "Annot") # This is the full name, with the 'organism' info

plots <- list()
tables <- list()
for (i in 1:length(fasta.file.names)) {
  fasta.file <- read.fasta(fasta.file.names[i], seqtype = "AA", set.attributes = T, whole.header = T)
  # Only some of the info is useful
  fasta <- lapply(fasta.file, function(x) {
    aa <- x[1]
    return(aa)
  })
  
  # Make a data frame
  df <- as.data.frame(unlist(fasta))
  df$Annot <- row.names(df)
  row.names(df) <- c(1:nrow(df))
  
  # Extract the species name for matching
  annot <- strsplit(df$Annot, "] ")
  annot2 <- lapply(annot, function(x) {
    species <- grep("organism", x)
    species <- x[species]
    return(species)
  })
  
  df$species <- annot2
  df$Species <- lapply(strsplit(unlist(df$species), "="), function(x) x[2])
  
  # Add diel info from sleepy_fish
  df$diel <- sleepy_fish$Diel_Pattern[match(df$Species, sleepy_fish$Species_name, nomatch = NA)]
  # df$diel2 <- resolved_names$diel2[match(df$Species, resolved_names$unique_name, nomatch = NA)]
  df$diel3 <- ifelse(df$diel == "diurnal", "diurnal", ifelse(df$diel == "nocturnal", "nocturnal", ifelse(df$diel == "crepuscular/diurnal", "diurnal", ifelse(df$diel == "crepuscular/nocturnal", "nocturnal", "unknown"))))
  
  
  # Make and plot the table
  cont.table <- table(df$`unlist(fasta)`, df$diel3)
  # mosaic(cont.table, shade = T)
  
  # maybe works better
  
  mosaicplot(cont.table, shade = T)
  cont.table
}

names(plots) <- fasta.file.names
names(tables) <- fasta.file.names








