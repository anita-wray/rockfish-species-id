### Microhaplot Tutorial ###

######################################
### INSTALL PACKAGES & LOAD FUNCTIONS

# packages_needed <- c("devtools", "rlang","ggplot2", "scales", "ggpubr", "microhaplot","dplyr", "remotes", 
#                      "data.table", "RcppCNPy", "tidyverse")
# 
# for(i in 1:length(packages_needed)){
#   if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
#   suppressWarnings(suppressMessages(library(packages_needed[i], character.only = TRUE)))
# }
# 
# if(!("microhaplot" %in% installed.packages())){
#   devtools::install_github("ngthomas/microhaplot", build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
# }
library(microhaplot)

#################################################################################
# provide a directory path to host the microhaplot Shiny app

WORKDIR <- "C:/Users/diana.baetscher/Work/git-repos/rockfish-species-id" ## CHANGE TO FIT YOUR NEEDS
METADIR <- "C:/Users/diana.baetscher/Work/git-repos/rockfish-species-id/data/rds_files"         ## CHANGE TO FIT YOUR NEEDS
PREFIX <- "sebastes"                                ## CHANGE TO FIT YOUR NEEDS
microhaplot::mvShinyHaplot(METADIR)
out.path<- tempdir()
app.path <- file.path(METADIR, "microhaplot")
microhaplot::runShinyHaplot(app.path)

#################################################################################
### Histograms


setwd(paste0(meta.path, "/idx_files/"))
ALLDATAFILENAMES <- Sys.glob("*")
allIDX_file_list <- as.list(ALLDATAFILENAMES)
idx_df <- allIDX_file_list %>%
  set_names(nm = ALLDATAFILENAMES) %>%
  map_dfr(
    ~ read_delim(.x, col_types = cols(), col_names = c("locus","length", "mapped", "unmapped"), delim = "\t")
  )
setwd(WORKDIR)


idx_grouped_df <- idx_df %>%
                      group_by(locus) %>%
                      summarize(sum(mapped)) %>%
                      rename(reads = `sum(mapped)`) %>%
                      arrange(desc(reads))

write.csv(idx_grouped_df, paste0(METADIR, PREFIX, "_baseline_part1_read_depth_per_locus.csv"), quote = FALSE, row.names = FALSE)


reads.perlocus_plot <- idx_grouped_df %>%
#  filter(reads) %>%
  ggplot() +
  geom_bar(aes(x = reorder(locus,reads), y = reads), stat = "identity") +
  ggtitle(paste0("Mean locus coverage = ", round(mean(idx_grouped_df$reads), 1))) +
  xlab("Locus") +
  theme(axis.text.x.bottom = element_blank())

jpeg(paste0(METADIR, PREFIX, "_baseline_part1_reads_perLocus.jpeg"), width = 15, height = 7, res = 150, units = "in")
plot(reads.perlocus_plot)
dev.off()

