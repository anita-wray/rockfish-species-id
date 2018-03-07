library(tidyverse)
library(readxl)
library(stringr)
library(lubridate)


source("R/rockfish-funcs2.R")


#### Call genos from the microhaplot RDS files ####

# this next step uses a little data frame that has the names and GTSeq run numbers
# I got that file like this:
# dsb:rds_files dianabaetscher$ pwd
# /Users/dianabaetscher/Desktop/NOAA_grad/git-repos/rockfish-species-id/data/rds_files
# dsb:rds_files dianabaetscher$ ls -l | awk 'BEGIN {print "gtseq_run", "file"} NR >1 && !/posinfo/ {num = $NF; gsub(/gtseq/, "", num); gsub(/[._].*$/, "", num);  print num, $NF}' > ../../data/rds-file-list.txt 

# get the names of the files
fdf <- read.table("data/rds-file-list.txt", stringsAsFactors = FALSE, header = TRUE) %>%
  tbl_df()
dir <- "data/rds_files"


# cycle over them, read them and add the gtseq_run column on each.
# at the end, bind them together.
genos_long <- lapply(1:nrow(fdf), function(i) {
  message("Working on ", fdf$file[i])
  call_genos_from_haplotRDS(path = file.path(dir, fdf$file[i])) %>%
    mutate(gtseq_run = fdf$gtseq_run[i]) %>%
    select(gtseq_run, everything())
}) %>%
  bind_rows()


# need to change the character GTseq to gtseq to match the sample sheet tibble
v1 <- genos_long %>%
  mutate(gtseq_run = replace(gtseq_run, str_detect(gtseq_run, "GTseq58"), "gtseq58"))
v2 <- v1 %>%
  mutate(gtseq_run = replace(gtseq_run, str_detect(gtseq_run, "GTseq59"), "gtseq59"))
v3 <- v2 %>%
  mutate(gtseq_run = replace(gtseq_run, str_detect(gtseq_run, "GTseq60"), "gtseq60"))
v4 <- v3 %>%
  mutate(gtseq_run = replace(gtseq_run, str_detect(gtseq_run, "GTseq62"), "gtseq62"))
genos_long <- v4 %>%
  mutate(gtseq_run = replace(gtseq_run, str_detect(gtseq_run, "GTseq63"), "gtseq63"))



# we go ahead and save it in data/processed, with xz compression
saveRDS(genos_long, file = "data/processed/called_genos.rds", compress = "xz")


#### Read the sample sheets from the Excel files ####

# Same drill here.  First we make a file that holds a data frame of file names:
# dsb:sample_sheets dianabaetscher$ pwd
# /Users/dianabaetscher/Desktop/NOAA_grad/git-repos/rockfish-species-id/data/sample_sheets
# dsb:sample_sheets dianabaetscher$ ls -l | awk 'BEGIN {print "gtseq_run", "file"} NR > 1 {num = $NF; gsub(/GTseq/, "", num); gsub(/_.*$/, "", num);  print num, $NF}' > ../../data/sample-sheet-file-list.txt

# There are a bunch of issues with GTseq62 - Hayley renamed the files
# and added additional samples to that directory.
# plus, apparently Lorne's samples don't have a NMFS DNA ID?

fdf <- read.table("data/sample-sheet-file-list.txt", stringsAsFactors = FALSE, header = TRUE) %>%
  tbl_df()
dir <- "data/sample_sheets"

sample_sheets <- lapply(1:nrow(fdf), function(i) {
  message("Working on ", fdf$file[i])
  read_excel(path = file.path(dir, fdf$file[i]), sheet = "sample_sheet", skip = 21) %>%
    filter(!is.na(Sample_Plate)) %>%  # read_excel sometimes reads a lot of blank lines...
    tidyr::separate(Sample_Plate, into = c("NMFS_DNA_ID", "ssBOX_ID", "ssBOX_POSITION")) %>%
    mutate(id = str_replace(Sample_ID, "s_0*", "s")) %>%
    mutate(gtseq_run = fdf$gtseq_run[i]) %>%
    select(gtseq_run, id, everything())
}) %>%
  bind_rows()

# save that.

saveRDS(sample_sheets, "data/processed/sample-sheet-tibble.rds", compress = "xz")



#### And finally, let's get the meta-data read in and cleaned up (if it needs it) ####
meta1 <- read_csv("data/nsf-rockfish-metadata.csv") %>%
  select(-`Marine::NMFS_DNA_ID`) %>%
  mutate(BATCH_ID = as.character(BATCH_ID),
         WEIGHT = as.numeric(WEIGHT)) # get rid of the duplicated ID column

# when we read that in, we lose the "None."s in the LEFTOVER_SAMPLE fields.  That 
# is OK for now.  

# that was the majority of the meta data, but there are two other sources,
# as we have some old SMURF data (plus a few other rockfish that are in that
# metadata set) and also fish from the multiple-paternity study with Sue Sogard,
# so we are going to get those meta data and bind them on there as well.
more_meta <- read_delim("data/more-rockfish-metadata.txt", delim = "\t") %>%
  select(-HATCHERY_MARK, -LANDFALL_PORT, -CRUISE, -HAUL)  %>%  # blow away a few of these columns that we don't need 
                                                          # so that the columns are congruent with meta
  select(-NMFS_DNA_ID_1)  # remove that duplicated ID column here too

names(meta1) <- str_replace(names(meta1), "^Marine::", "")  # get the "Marine::" out of some of the column names

meta <- bind_rows(meta1, more_meta) %>%
  mutate(COLLECTION_DATE = mdy(COLLECTION_DATE),
         PICK_DATE = mdy(PICK_DATE))

saveRDS(meta, "data/processed/meta-data-tibble.rds", compress = "xz")

# for re-doing this analysis, I just read in the existing sample sheet tibble and meta data tibble
#meta <- readRDS("data/processed/meta-data-tibble.rds")
#sample_sheets <- readRDS("data/processed/sample-sheet-tibble.rds")


#### In the end, let us get a data frame that includes genotypes for all the individuals  ####
# and which explicity has NAs in places where data are missing, and also 
# has the NMFS_DNA_ID on there
genos_long_explicit_NAs <- sample_sheets %>%
  select(gtseq_run, id, NMFS_DNA_ID) %>%
  unite(col = gid, sep = "_", gtseq_run, id, NMFS_DNA_ID) %>%
  select(gid) %>%
  unlist() %>%
  unname() %>%
  expand.grid(gid = ., locus = unique(genos_long$locus), gene_copy = 1:2, stringsAsFactors = FALSE) %>%
  tbl_df() %>% 
  separate(gid, into = c("gtseq_run", "id", "NMFS_DNA_ID"), convert = TRUE) %>%
  left_join(., genos_long) %>%
  arrange(gtseq_run, id, locus, gene_copy)

genos_long_explicit_NAs


# and then save that
saveRDS(genos_long_explicit_NAs, file = "data/processed/called_genos_na_explicit.rds", compress = "xz")
  
