unknown-species-id-template
================

28 February 2023

Walkthrough of species id for unknown samples:

Using the rds file output from microhaplot: 1. read in rds files 2.
apply read depth filters 3. apply allele balance filter

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.4.1      v purrr   0.3.4 
    ## v tibble  3.1.2      v dplyr   1.0.10
    ## v tidyr   1.2.0      v stringr 1.4.0 
    ## v readr   1.4.0      v forcats 0.5.1

    ## Warning: package 'ggplot2' was built under R version 4.1.3

    ## Warning: package 'tidyr' was built under R version 4.1.3

    ## Warning: package 'dplyr' was built under R version 4.1.3

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(readxl)
library(stringr)
library(lubridate)
```

    ## 
    ## Attaching package: 'lubridate'

    ## The following objects are masked from 'package:base':
    ## 
    ##     date, intersect, setdiff, union

``` r
library(rubias)
```

    ## Warning: package 'rubias' was built under R version 4.1.3

``` r
library(ggpattern)
```

    ## Warning: package 'ggpattern' was built under R version 4.1.3

``` r
source("../R/rockfish-funcs2.R")


#### Call genos from the microhaplot rds files ####

# this next step uses a little data frame that has the names and GTSeq run numbers
# I got that file like this:
# dsb:rds_files dianabaetscher$ pwd
# /Users/dianabaetscher/Desktop/NOAA_grad/git-repos/rockfish-species-id/new_baseline_data/feather_files
# dsb:feather_files dianabaetscher$ ls -l | awk 'BEGIN {print "gtseq_run", "file"} NR >1 && !/posinfo/ {num = $NF; gsub(/gtseq/, "", num); gsub(/[._].*$/, "", num);  print num, $NF}' > ../../data/rds-file-list.txt 

# get the names of the files
fdf <- read.table("../data_AFSC/rds_files/rds-file-list.txt", stringsAsFactors = FALSE, header = TRUE) %>%
  tbl_df()
```

    ## Warning: `tbl_df()` was deprecated in dplyr 1.0.0.

    ## Warning: Please use `tibble::as_tibble()` instead.

``` r
dir <- "../data_AFSC/rds_files/select_for_analysis/"


# cycle over them, read them and add the gtseq_run column on each.
# at the end, bind them together.
genos_long <- lapply(1:nrow(fdf), function(i) {
  message("Working on ", fdf$file[i])
  call_genos_from_haplotRDS(path = file.path(dir, fdf$file[i])) %>%
    mutate(gtseq_run = fdf$gtseq_run[i]) %>%
    select(gtseq_run, everything())
}) %>%
  bind_rows()
```

    ## Working on sebastes_GOA_larvae.rds

    ## Joining, by = c("id", "locus", "rank")Working on sebastes_gtseq_test1.rds
    ## Joining, by = c("id", "locus", "rank")Working on sebastes_mixedPlates.rds
    ## Joining, by = c("id", "locus", "rank")Working on sebastes_GOA_larvae2.rds
    ## Joining, by = c("id", "locus", "rank")

``` r
# we go ahead and save it in data/processed, with xz compression
saveRDS(genos_long, file = "../data_AFSC/processed/called_genos.rds", compress = "xz")


#### In the end, let us get a data frame that includes genotypes for all the individuals  ####
# and which explicitly has NAs in places where data are missing, and also 
# has the NMFS_DNA_ID on there
genos_long_explicit_NAs <- genos_long %>%
  select(gtseq_run, id) %>%
  unique() %>%
  unite(col = gid, sep = "_", gtseq_run, id) %>%
  select(gid) %>%
  unlist() %>%
  unname() %>%
  expand.grid(gid = ., locus = unique(genos_long$locus), gene_copy = 1:2, stringsAsFactors = FALSE) %>%
  tbl_df() %>% 
  separate(gid, into = c("gtseq_run", "id"), convert = TRUE) %>%
  left_join(., genos_long) %>%
  arrange(gtseq_run, id, locus, gene_copy)
```

    ## Joining, by = c("gtseq_run", "id", "locus", "gene_copy")

``` r
# and then save that
saveRDS(genos_long_explicit_NAs, file = "../data_AFSC/processed/called_genos_na_explicit.rds", compress = "xz")
```

Using those genotypes…

``` r
genos_long_explicit_NAs %>%
  group_by(gtseq_run, id) %>%
  tally() 
```

    ## # A tibble: 568 x 3
    ## # Groups:   gtseq_run [4]
    ##    gtseq_run id        n
    ##    <chr>     <chr> <int>
    ##  1 gtseq1    s100    180
    ##  2 gtseq1    s101    180
    ##  3 gtseq1    s102    180
    ##  4 gtseq1    s103    180
    ##  5 gtseq1    s104    180
    ##  6 gtseq1    s105    180
    ##  7 gtseq1    s106    180
    ##  8 gtseq1    s107    180
    ##  9 gtseq1    s108    180
    ## 10 gtseq1    s109    180
    ## # ... with 558 more rows

381 samples

Look at missing data: 180 gene copies total (90 loci x2)

``` r
ind_to_toss <- genos_long_explicit_NAs %>%
  group_by(gtseq_run, id) %>%
  filter(is.na(allele)) %>% # missing data
  tally() %>%
  arrange(desc(n)) %>% # remove samples with >20% missing data
  filter(n > 36)

# remove those from the df
genos_ind_filtered <- genos_long_explicit_NAs %>%
  anti_join(., ind_to_toss)
```

    ## Joining, by = c("gtseq_run", "id")

Load baseline data

``` r
# baseline data - curated, 997 indivs
baseline <- readRDS("../new_baseline_data/processed/sebastes_spp_id_baseline_haplotypes.rds")

# remove the 6 loci that had HWE and other issues
to_remove <- read_csv("../data/loci_to_remove.csv")
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   locus = col_character()
    ## )

``` r
baseline90 <- baseline %>%
  anti_join(., to_remove)
```

    ## Joining, by = "locus"

``` r
# remind myself which species are in the baseline:
baseline90 %>%
  select(collection) %>%
  unique() %>%
  arrange()
```

    ## # A tibble: 54 x 1
    ##    collection   
    ##    <chr>        
    ##  1 aleutianus   
    ##  2 alutus       
    ##  3 auriculatus  
    ##  4 aurora       
    ##  5 babcocki     
    ##  6 borealis     
    ##  7 caurinus     
    ##  8 chlorostictus
    ##  9 constellatus 
    ## 10 crameri      
    ## # ... with 44 more rows

``` r
tossers <- baseline90 %>%
  select(indiv, gtseq_run, id) %>%
  unique() %>%
  group_by(indiv) %>%
  tally() %>%
  filter(n >1)

baseline90_one_each <- baseline90 %>%
  anti_join(., tossers)
```

    ## Joining, by = "indiv"

``` r
# baseline data - curated, 997 indivs
baseline_spp_info <- baseline90_one_each %>%
  select(sample_type, repunit, collection, indiv, gtseq_run, id, species) %>%
  unique()
baseline_spp_info$gtseq_run <- as.character(baseline_spp_info$gtseq_run)
```

``` r
# slim that down to just the matching field with the unknowns
for_alleidx <- baseline90_one_each %>%
  select(-indiv, -c(1:3, 12:13), -species)
  

for_alleidx$gtseq_run <- as.character(for_alleidx$gtseq_run)
```

``` r
# merge the two dataframes
merged_df <- bind_rows(for_alleidx, genos_ind_filtered)

# first make integers of the alleles
alle_idxs <- merged_df %>% 
  dplyr::select(gtseq_run, id, locus, gene_copy, allele) %>%
  group_by(locus) %>%
  mutate(alleidx = as.integer(factor(allele, levels = unique(allele)))) %>%
  ungroup() %>%
  arrange(gtseq_run, id, locus, alleidx) # rubias can handle NA's, so no need to change them to 0's

  
# and spread the alleles
two_col <- alle_idxs %>%
  #group_by(indiv, locus) %>%
  unite(loc, locus, gene_copy, sep = ".") %>%
  #ungroup() %>%
  select(-allele) %>%
  pivot_wider(names_from = loc, values_from = alleidx) 
```

add back on info for reference and make two-column format for rubias

``` r
# baseline
reference <- two_col %>%
  left_join(., baseline_spp_info) %>%
  filter(!is.na(species)) %>%
  select(-gtseq_run, -id, -species) %>%
  select(sample_type, repunit, collection, indiv, everything())
```

    ## Joining, by = c("gtseq_run", "id")

``` r
# mixture
rubias_mix <- two_col %>%
  anti_join(., baseline_spp_info) %>%
  mutate(sample_type = "mixture", collection = "larvae", repunit = NA) %>%
  select(sample_type, repunit, collection, everything()) %>%
  unite(gtseq_run, id, col = "indiv", sep = "_")
```

    ## Joining, by = c("gtseq_run", "id")

## Mixture assignment with rubias

``` r
rubias_output <- infer_mixture(reference = reference, mixture = rubias_mix, gen_start_col = 5)
```

    ## Collating data; compiling reference allele frequencies, etc.   time: 1.11 seconds
    ## Computing reference locus specific means and variances for computing mixture z-scores   time: 0.14 seconds
    ## Working on mixture collection: larvae with 409 individuals
    ##   calculating log-likelihoods of the mixture individuals.   time: 0.13 seconds
    ##   performing 2000 total sweeps, 100 of which are burn-in and will not be used in computing averages in method "MCMC"   time: 0.35 seconds
    ##   tidying output into a tibble.   time: 0.14 seconds

``` r
# take the top output for each sample
top_assign <- rubias_output$indiv_posteriors %>%
  group_by(indiv) %>%
  slice_max(., order_by = PofZ)


top_assign
```

    ## # A tibble: 409 x 10
    ## # Groups:   indiv [409]
    ##    mixture~1 indiv repunit colle~2  PofZ log_l~3 z_score n_non~4 n_mis~5 missi~6
    ##    <chr>     <chr> <chr>   <chr>   <dbl>   <dbl>   <dbl>   <int>   <int> <list> 
    ##  1 larvae    gtse~ rufina~ rufina~     1   -464.   -57.7      72      18 <int>  
    ##  2 larvae    gtse~ rufina~ rufina~     1   -477.   -59.0      74      16 <int>  
    ##  3 larvae    gtse~ rufina~ rufina~     1   -463.   -56.1      73      17 <int>  
    ##  4 larvae    gtse~ rufina~ rufina~     1   -494.   -57.7      77      13 <int>  
    ##  5 larvae    gtse~ rufina~ rufina~     1   -486.   -59.7      76      14 <int>  
    ##  6 larvae    gtse~ rufina~ rufina~     1   -487.   -63.9      76      14 <int>  
    ##  7 larvae    gtse~ rufina~ rufina~     1   -486.   -59.7      76      14 <int>  
    ##  8 larvae    gtse~ rufina~ rufina~     1   -448.   -55.5      72      18 <int>  
    ##  9 larvae    gtse~ rufina~ rufina~     1   -464.   -57.0      73      17 <int>  
    ## 10 larvae    gtse~ rufina~ rufina~     1   -476.   -55.6      74      16 <int>  
    ## # ... with 399 more rows, and abbreviated variable names 1: mixture_collection,
    ## #   2: collection, 3: log_likelihood, 4: n_non_miss_loci, 5: n_miss_loci,
    ## #   6: missing_loci

Check on z-scores:

``` r
top_assign %>%
  ggplot(aes(x = z_score)) +
  geom_histogram()
```

    ## 

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](05-unknown-species-id-template_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggsave("pdf_outputs/larvae_mixed_unknowns_Zscores.pdf")
```

    ## 
    ## Saving 7 x 5 in image
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

These z-scores suggest that the larvae are not only species not in the
baseline, but some of them are very phylogenetically different from the
baseline species.

Some of that assessment depends on what the cherrypicked plates look
like:

``` r
top_assign %>%
  filter(str_detect(indiv, "gtseq2")) %>%
  filter(PofZ > 0.95) %>% # none of the -60 z-scores in this batch
  ggplot(aes(x = z_score)) +
  geom_histogram()
```

    ## 

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](05-unknown-species-id-template_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Read in some metadata to take a closer look at those
samples/assignments.

``` r
POPmix_samplesheet <- read_csv("../data_AFSC/samplesheets/20230110_POPMix_POP_PCOD_BFAL.csv", skip = 19)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   Sample_ID = col_character(),
    ##   Sample_Plate = col_character(),
    ##   Sample_Well = col_character(),
    ##   I7_Index_ID = col_character(),
    ##   index = col_character(),
    ##   I5_Index_ID = col_character(),
    ##   index2 = col_character(),
    ##   Sample_Project = col_character(),
    ##   Description = col_logical()
    ## )

``` r
meta_spp <- read_csv("../data_AFSC/POP_northerns_duskies_to_gtseq_FIXEDmetadata20230327.csv")
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   .default = col_character(),
    ##   ABLG = col_double(),
    ##   CaptureTime = col_logical(),
    ##   Comments = col_logical(),
    ##   DepthM = col_double(),
    ##   DNAngul = col_double(),
    ##   ExtractionPlate = col_double(),
    ##   FishDepthMM = col_logical(),
    ##   LengthMM = col_double(),
    ##   ReleaseTime = col_logical(),
    ##   StartLatitudeDD = col_double(),
    ##   StartLongitudeDD = col_double(),
    ##   TempC = col_double(),
    ##   WeightG = col_double()
    ## )
    ## i Use `spec()` for the full column specifications.

``` r
meta_spp$ABLG <- as.character(meta_spp$ABLG)

# join those to get the species
meta_for_gtseqMix <- POPmix_samplesheet %>%
  filter(Sample_Project == "POPMix") %>%
       left_join(., meta_spp, by = c("Sample_ID" = "ABLG")) %>%
  select(Sample_ID, SpeciesName, AlternateID_s_, Sample_Plate, Sample_Well) %>%
  mutate(id = row_number()) %>%
  mutate(s = "s") %>%
  unite(col = "s_num", s, id, sep = "")
```

``` r
# POP and duskies mixed plate
spp_id_w_knowns <- top_assign %>%
  filter(str_detect(indiv, "gtseq2")) %>%
  separate(indiv, into = c("gtseq_run", "id")) %>%
  left_join(., meta_for_gtseqMix, by = c("id" = "s_num")) %>%
  select(collection, SpeciesName, z_score, everything()) %>%
  rename(rubias_assignment = collection)

spp_id_w_knowns %>%
  ggplot(aes(x = z_score, fill = SpeciesName)) +
  geom_histogram(alpha = 0.3) +
  theme_bw()
```

    ## 

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](05-unknown-species-id-template_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
spp_id_w_knowns %>%
  ggplot(aes(x = z_score, fill = SpeciesName)) +
  geom_histogram(alpha = 0.5) +
  facet_grid(rows = vars(SpeciesName)) +
  theme_bw() +
  theme(
    strip.text.y = element_text(size = 8)) +
  labs(title = "POP and dusky mixed plate - species ID")
```

    ## 
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](05-unknown-species-id-template_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
ggsave("pdf_outputs/known_sebastes_POP_dusky_northern_sppID.pdf", height = 5, width = 5)
```

    ## 
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

``` r
spp_id_w_knowns %>%
  select(-missing_loci) %>%
  write_csv("csv_outputs/known_northern_POP_dusky_rubias_id.csv")

spp_id_w_knowns %>%
  filter(SpeciesName == "Sebastes polyspinis")
```

    ## # A tibble: 48 x 16
    ##    rubias_~1 Speci~2 z_score mixtu~3 gtseq~4 id    repunit  PofZ log_l~5 n_non~6
    ##    <chr>     <chr>     <dbl> <chr>   <chr>   <chr> <chr>   <dbl>   <dbl>   <int>
    ##  1 polyspin~ Sebast~ -0.108  larvae  gtseq2  s1    polysp~     1   -40.2      82
    ##  2 polyspin~ Sebast~  0.816  larvae  gtseq2  s10   polysp~     1   -33.0      80
    ##  3 polyspin~ Sebast~  0.0213 larvae  gtseq2  s11   polysp~     1   -39.3      82
    ##  4 polyspin~ Sebast~ -2.95   larvae  gtseq2  s12   polysp~     1   -60.6      84
    ##  5 polyspin~ Sebast~  0.976  larvae  gtseq2  s13   polysp~     1   -33.3      83
    ##  6 polyspin~ Sebast~  0.712  larvae  gtseq2  s14   polysp~     1   -35.2      83
    ##  7 polyspin~ Sebast~ -1.98   larvae  gtseq2  s15   polysp~     1   -53.8      84
    ##  8 polyspin~ Sebast~ -0.835  larvae  gtseq2  s16   polysp~     1   -45.2      81
    ##  9 polyspin~ Sebast~ -1.55   larvae  gtseq2  s17   polysp~     1   -50.9      83
    ## 10 polyspin~ Sebast~ -0.576  larvae  gtseq2  s18   polysp~     1   -44.1      83
    ## # ... with 38 more rows, 6 more variables: n_miss_loci <int>,
    ## #   missing_loci <list>, Sample_ID <chr>, AlternateID_s_ <chr>,
    ## #   Sample_Plate <chr>, Sample_Well <chr>, and abbreviated variable names
    ## #   1: rubias_assignment, 2: SpeciesName, 3: mixture_collection, 4: gtseq_run,
    ## #   5: log_likelihood, 6: n_non_miss_loci

That actually makes me feel considerably better about the assignment.
The z-scores for Northerns look normal (and Northerns are panmictic)
whereas the POP has a wider distribution (because of pop structure,
likely) and the duskies are the first set of z-score outliers.

In fact, we could use the z-scores to some extent, to give us a hint
about the larval samples.

``` r
goa <- top_assign %>%
  filter(!str_detect(indiv, "gtseq2")) 

# combined plot
ggplot() +
  geom_histogram(data = goa, aes(x = z_score), alpha = 0.5) +
  geom_histogram(data = spp_id_w_knowns, aes(x = z_score, fill = SpeciesName), alpha = 0.25) +
  theme_bw() +
  labs(title = "GOA unknown larvae - Sebastes species ID")
```

    ## 

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](05-unknown-species-id-template_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
ggsave("pdf_outputs/GOA_unknowns_w_unknownSpp2.pdf")
```

    ## 
    ## Saving 7 x 5 in image
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

What was the species assignment for the -60 bunch?

``` r
goa %>%
  filter(z_score < -20)
```

    ## # A tibble: 13 x 10
    ## # Groups:   indiv [13]
    ##    mixture~1 indiv repunit colle~2  PofZ log_l~3 z_score n_non~4 n_mis~5 missi~6
    ##    <chr>     <chr> <chr>   <chr>   <dbl>   <dbl>   <dbl>   <int>   <int> <list> 
    ##  1 larvae    gtse~ rufina~ rufina~     1   -464.   -57.7      72      18 <int>  
    ##  2 larvae    gtse~ rufina~ rufina~     1   -477.   -59.0      74      16 <int>  
    ##  3 larvae    gtse~ rufina~ rufina~     1   -463.   -56.1      73      17 <int>  
    ##  4 larvae    gtse~ rufina~ rufina~     1   -494.   -57.7      77      13 <int>  
    ##  5 larvae    gtse~ rufina~ rufina~     1   -486.   -59.7      76      14 <int>  
    ##  6 larvae    gtse~ rufina~ rufina~     1   -487.   -63.9      76      14 <int>  
    ##  7 larvae    gtse~ rufina~ rufina~     1   -486.   -59.7      76      14 <int>  
    ##  8 larvae    gtse~ rufina~ rufina~     1   -448.   -55.5      72      18 <int>  
    ##  9 larvae    gtse~ rufina~ rufina~     1   -464.   -57.0      73      17 <int>  
    ## 10 larvae    gtse~ rufina~ rufina~     1   -476.   -55.6      74      16 <int>  
    ## 11 larvae    gtse~ rufina~ rufina~     1   -455.   -55.0      72      18 <int>  
    ## 12 larvae    gtse~ rufina~ rufina~     1   -473.   -59.0      72      18 <int>  
    ## 13 larvae    gtse~ rufina~ rufina~     1   -489.   -60.1      77      13 <int>  
    ## # ... with abbreviated variable names 1: mixture_collection, 2: collection,
    ## #   3: log_likelihood, 4: n_non_miss_loci, 5: n_miss_loci, 6: missing_loci

``` r
# those samples with normal parameters
normal_Zs <- top_assign %>%
  filter(z_score > -3 & z_score < 3) 

normal_Zs %>%
  ggplot(aes(x = z_score)) +
  geom_histogram()
```

    ## 

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](05-unknown-species-id-template_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

How many larvae from which species in the “normal” range for baseline
species?

``` r
normal_Zs %>%
  group_by(collection) %>%
  tally()
```

    ## # A tibble: 5 x 2
    ##   collection       n
    ##   <chr>        <int>
    ## 1 alutus          89
    ## 2 babcocki         1
    ## 3 maliger         11
    ## 4 nigrocinctus     2
    ## 5 polyspinis      49

POP, northerns, quillback, redbanded, and tiger rockfishes. All
reasonable.

## Check self-assignment of baseline

Sanity-check

Look at self-assignment:

``` r
self_test <- self_assign(reference = reference, gen_start_col = 5)
```

    ## Summary Statistics:
    ## 
    ## 992 Individuals in Sample
    ## 
    ## 90 Loci: Plate_1_A01_Sat_GW603857_consensus.1, Plate_1_A11_Sat_GE820299_consensus.1, Plate_2_A09_Sat_EW986980_consensus.1, Plate_2_C08_Sat_EW987116_consensus.1, Plate_2_G06_Sat_EW987118_consensus.1, Plate_3_C03_Sat_GE798118_consensus.1, Plate_4_E10_Sat_EW976030_consensus.1, Plate_4_G06_Sat_EW976181_consensus.1, tag_id_1049.1, tag_id_108.2, tag_id_1184.1, tag_id_1229.1, tag_id_1272.1, tag_id_1366.1, tag_id_1428.1, tag_id_143.1, tag_id_1441.1, tag_id_1449.1, tag_id_1471.1, tag_id_1498.1, tag_id_1558.1, tag_id_1576.1, tag_id_1598.1, tag_id_1604.1, tag_id_1613.1, tag_id_162.1, tag_id_1652.1, tag_id_170.1, tag_id_1708.1, tag_id_1748.1, tag_id_1751.1, tag_id_1762.1, tag_id_179.1, tag_id_1804.1, tag_id_1808.1, tag_id_1810.1, tag_id_1836.1, tag_id_1850.1, tag_id_1880.1, tag_id_1889.1, tag_id_1915.1, tag_id_1950.1, tag_id_1961.1, tag_id_1966.2, tag_id_1982.1, tag_id_1994.1, tag_id_1999.1, tag_id_2008.1, tag_id_2009.2, tag_id_2017.1, tag_id_2062.1, tag_id_2082.1, tag_id_2114.1, tag_id_2134.1, tag_id_2155.1, tag_id_2178.1, tag_id_2182.1, tag_id_220.1, tag_id_2203.1, tag_id_221.1, tag_id_2214.1, tag_id_2237.1, tag_id_2247.1, tag_id_2258.1, tag_id_2301.1, tag_id_2319.1, tag_id_2368.1, tag_id_2499.1, tag_id_250.1, tag_id_2607.1, tag_id_2635.1, tag_id_265.1, tag_id_325.1, tag_id_402.1, tag_id_410.1, tag_id_436.1, tag_id_55.1, tag_id_572.1, tag_id_67.1, tag_id_770.1, tag_id_788.1, tag_id_843.1, tag_id_855.1, tag_id_874.1, tag_id_875.1, tag_id_879.1, tag_id_913.1, tag_id_942.1, tag_id_981.1, tag_id_987.1
    ## 
    ## 53 Reporting Units: melanops, caurinus, hopkinsi, mystinus, atrovirens, chrysomelas, auriculatus, entomelas, jordani, maliger, miniatus, nebulosus, goodei, flavidus, chlorostictus, diaconus, elongatus, ovalis, paucispinis, pinniger, rastrelliger, saxicola, serranoides, crameri, proriger, rosaceus, wilsoni, diploproa, semicinctus, babcocki, reedi, dallii, aurora, umbrosus, levis, constellatus, oculatus, melanostomus, zacentrus, ruberrimus, rubrivinctus, serriceps, rufus, alutus, emphaeus, nigrocinctus, ensifer, moseri, rufinanus, melanostictus, aleutianus, borealis, polyspinis
    ## 
    ## 54 Collections: melanops, caurinus, hopkinsi, mystinus, atrovirens, chrysomelas, carnatus, auriculatus, entomelas, jordani, maliger, miniatus, nebulosus, goodei, flavidus, chlorostictus, diaconus, elongatus, ovalis, paucispinis, pinniger, rastrelliger, saxicola, serranoides, crameri, proriger, rosaceus, wilsoni, diploproa, semicinctus, babcocki, reedi, dallii, aurora, umbrosus, levis, constellatus, oculatus, melanostomus, zacentrus, ruberrimus, rubrivinctus, serriceps, rufus, alutus, emphaeus, nigrocinctus, ensifer, moseri, rufinanus, melanostictus, aleutianus, borealis, polyspinis
    ## 
    ## 8.55% of allelic data identified as missing

``` r
self_test %>%
  group_by(indiv) %>%
  slice_max(., order_by = scaled_likelihood) %>%
  ggplot(aes(x = z_score)) +
  geom_histogram()
```

    ## 

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](05-unknown-species-id-template_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
ggsave("pdf_outputs/baseline_zscores.pdf")
```

    ## 
    ## Saving 7 x 5 in image
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

Self-assignment looks fine.

I’m beginning to think that what we’re seeing in the z-scores is (1) the
population structure in POP, (2) duskies missing from the baseline, and
(3) some weird species - maybe not sebastes?? in the larval collections.

##### DO NOT MODIFY BELOW THIS LINE

## New primer pool

Take a look at the distribution of reads across loci:

``` r
loc_depth <- genos_long_explicit_NAs %>%
  group_by(gtseq_run, id, locus) %>%
  filter(!is.na(depth)) %>% # remove the NAs, which interfere with sum
  ungroup() %>%
  group_by(locus) %>%
  summarise(mean_depth = mean(depth)) %>%
  arrange(desc(mean_depth)) 

loc_depth %>%
  ggplot(aes(x = reorder(locus, mean_depth), y = mean_depth)) +
  geom_bar(stat = "identity") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.95)
  )
```

![](05-unknown-species-id-template_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Try changing some primer concentrations:

C1\* V1 = C2 \*V2

``` r
primer_conc <- loc_depth %>%
  arrange(desc(mean_depth)) %>%
  mutate(primer_multiplier = ifelse(mean_depth < 100, 2, 1)) %>%
  mutate(primer_multiplier = ifelse(mean_depth > 500, 0.5, primer_multiplier)) %>%
  mutate(primer_conc = 0.25) %>%
  mutate(new_primer_conc = primer_multiplier*primer_conc) %>%
  mutate(C1 = 200) %>%
  mutate(V2 = 1500) %>%
  mutate(V1 = round(V2*new_primer_conc/C1)) %>%
  rename(vol_to_add = V1, total_vol = V2)
  # summarise(sum(vol_to_add))
```

202 ul of primers

``` r
1500-202 
```

    ## [1] 1298

1298 ul of Tris or water to add.

Read in primer order for plate map:

``` r
primerOrder <- read_xlsx("../../../GTseq/sebastesSppID/PrimerOrder.xlsx")

primer_conc$locus <- gsub("tag_id", "Sat", primer_conc$locus)

plate_primers <- primer_conc %>%
  filter(str_detect(locus, "Plate_")) %>%
  separate(locus, into = c("plate", "n", "well", "spp", "label", "drop")) %>%
  unite(col = "locus", spp, label) %>%
  select(-c(1:3, 5))



primerOrder$`Sequence Name` <- gsub("_F", "", primerOrder$`Sequence Name`)
order_short <- primerOrder %>%
  select(`Well Position`, `Sequence Name`)


primer_conc %>%
  filter(!str_detect(locus, "Plate_")) %>%
  bind_rows(., plate_primers) %>%
  left_join(., order_short, by = c("locus" = "Sequence Name")) %>%
  unique() %>%
   write_csv("csv_outputs/sebastes_spp_id_primerPool_20230228.csv")
```
