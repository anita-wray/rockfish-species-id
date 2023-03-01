species identification in rockfishes
================
19 December, 2022

-   <a href="#preliminaries" id="toc-preliminaries">Preliminaries</a>
-   <a href="#data-files" id="toc-data-files">Data files</a>
-   <a href="#data-filtering" id="toc-data-filtering">Data filtering</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

This is the repository with data and analyses for rockfish species
identification. This project arose from ambiguous identifications of
larval rockfishes in Monterey Bay, California, with expanded coverage of
southern and more northern species thanks to collaborators within NOAA
science centers.

We use amplicon-sequencing on our MiSeq instrument to generate genotype
data for samples. Using these samples, we created a vcf file that
contains variant sites from 54 species of rockfishes, mostly from the
Northeast Pacific.

`sebastes_sppID_combined_filtered.recode.vcf`

## Preliminaries

In order to perform Sebastes species ID…

1.  Process the GTseq run with the samples you wish to identify using
    the combined Sebastes vcf in MICROHAPLOT. Scripts for processing raw
    data through microhaplot are available here:

<https://github.com/AFSC-Genetics/GTseq_microhaplot>

2.  Open .rds files in the microhaplot shiny app. Download diploid
    genotype tables to import into R.

## Data files

1.  .csv files created by the R software program MICROHAPLOT should be
    placed in the `Rmd_AFSC/microhaplot_outputs` directory, these data
    are read into R in `02-test-pca-w-unknowns.Rmd`.

2.  The MiSeq sample sheets associated with these genotype data should
    be placed in `data/sample_sheets`. This is the means by which sample
    ID numbers from the sequencing runs can be connected back to NMFS
    DNA ID as recorded in meta data and the repository.

## Data filtering

1.  Remove the 6 loci from the baseline that were most prone to missing
    data and/or often have more than 3 haplotypes passing filter across
    individuals:

-   tag_id_1166 (often \>2 haplotypes)

-   tag_id_934 (missing data and sometimes \>2 haplotypes)

-   tag_id_2513 (appears to be some sort of repetitive element for some
    species outside of KGBC)

-   tag_id_1871 (failed in many individuals)

-   tag_id_1399 (failed in many individuals)

-   tag_id_914 (failed in many individuals)

2.  Typically we institute a missing data threshold of \~10% for
    excluding samples, but such a criterion doesn’t make sense given the
    dramatic ascertainment bias from kelp rockfish to less closely
    related species across phylogenetic distance. We can explore the
    missing data criterion in more depth, since there might be fixed
    alleles that make even a small number of loci valid for species ID.
    
## Species assignment
1.  An example of species identification of a mixture of samples using
    rubias is outlined in `Rmd_AFSC/05-unknown-species-id-template.Rmd`.
    
2. High and low z-scores are indicative of species (or populations of species) that were not included in the reference baseline.

## Updating the baseline
1. The reference baseline can be updated by genotyping additional species with the same set of primers, using the same bioinformatic workflow, and then merging VCF files from the prior baseline VCF.

2. The newly merged VCF file can then be called for analysis in microhaplot.
