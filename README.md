species identification in rockfishes
================
11 May, 2018

-   [Preliminaries](#preliminaries)
-   [Data files](#data-files)
-   [Data filtering](#data-filtering)

<!-- README.md is generated from README.Rmd. Please edit that file -->
This is an Rstudio project with data and analyses for species identification work that arose from the NSF rockfish dispersal project. Hayley helped expand the number of species sequenced, and I am creating a more reproducible and streamlined process here.

We use amplicon-sequencing on our MiSeq instrument to generate genotype data for samples. Using these samples, we created a vcf file that contains variant sites from 49 species of rockfishes, mostly from the Northeast Pacific. These data correspond to GT-seq runs 11, 12, 19, 48, 54, 55, and 62 and is located in this directory on Megabox.

`/home/biopipe/Genetics_Lab_Data/GTseq/sebastes09282017/map/seb_SppID_all96loci_sorted_10022017.vcf`

I have modified Hayley's original workflow to adapt it to an R project and make it more reproducible: this is still a work in progress.

Preliminaries
-------------

In order to perform Sebastes species ID...

1.  Process the GTseq run with the samples you wish to identify using the combined Sebastes vcf in MICROHAPLOT. The most current vcf can be found at:

/home/biopipe/Genetics\_Lab\_Data/GTseq/sebastes09282017/map/seb\_SppID\_all96loci\_sorted\_10022017.vcf

1.  Copy rds files from HaPLOType (the minimum read depth and min allelic ratio filters are applied in `rockfish-funcs2.R` when reading the rds files into this R project)

Data files
----------

1.  rds files created by the R software program MICROHAPLOT should be placed in the `data/rds_files` directory, with instructions for reading in genotypes data in `R-main/01-compile-and-tidy-data.R`.

2.  The MiSeq sample sheets associated with these genotype data should be placed in `data/sample_sheets`. This is the means by which sample ID numbers from the sequencing runs can be connected back to NMFS DNA ID as recorded in meta data and the repository.

Data filtering
--------------

1.  Remove the 6 loci that were most prone to missing data and/or often have more than 3 haplotypes passing filter across individuals:

-   tag\_id\_1166 (often &gt;2 haplotypes)

-   tag\_id\_934 (missing data and sometimes &gt;2 haplotypes)

-   tag\_id\_2513 (appears to be some sort of repetitive element for some species outside of KGBC)

-   tag\_id\_1871 (failed in many individuals)

-   tag\_id\_1399 (failed in many individuals)

-   tag\_id\_914 (failed in many individuals)

1.  Typically we institute a missing data threshold of ~10% for excluding samples, but such a criterion doesn't make sense given the dramatic ascertainment bias from kelp rockfish to less closely related species across phylogenetic distance. For now, I am keeping all individuals in the baseline.

2.  An example of species identification of a mixture of samples using rubias is outlined in `02-assign-mixture-w-rubias.Rmd`.
