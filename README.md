species identification in rockfishes
================
23 February, 2018

-   [Hayley's workflow](#hayleys-workflow)

<!-- README.md is generated from README.Rmd. Please edit that file -->
This is an Rstudio project with data and analyses for species identification work that arose from the NSF rockfish dispersal project. Hayley helped expand the number of species sequenced, and I am creating a more reproducible and streamlined process here.

We use amplicon-sequencing on our MiSeq instrument to generate genotype data for samples. Using these samples, we created a vcf file that contains variant sites from 49 species of rockfishes, mostly from the Northeast Pacific. These data correspond to GT-seq runs 11, 12, 19, 48, 54, 55, and 62 (maybe - need to check this with Hayley) and is located in this directory on Megabox. `/home/biopipe/Genetics_Lab_Data/GTseq/sebastes09282017/map/seb_SppID_all96loci_sorted_10022017.vcf`

For now, this is how the species identification has been done...

Hayley's workflow
-----------------

Workflow for performing Sebastes Species ID

1.  Run the GTseq run with the samples you wish to ID against the combined Sebastes vcf in HaPLOT. The most current vcf can be found at: /home/biopipe/Genetics\_Lab\_Data/GTseq/sebastes09282017/map/seb\_SppID\_all96loci\_sorted\_10022017.vcf

2.  Visualize in HaPLOType and download reported individual haplotype file after applying your preferred filtration criteria (i.e. min read depth = 10, min allelic ratio = 0.2)

3.  Combine this haplotype report file to the haplotype report that includes all individuals in the baseline (i.e. all of the individuals that contributed to the vcf). The file with all individuals in the baseline is “reported\_haplotype\_SebSppID\_11102017.csv” You may need to change the sampleID of the individuals on your GTseq run so that they do not overlap with samples in the baseline. Make sure to keep track of how samples were re-numbered.

4.  Create new R project for this run and use script to convert your GTseq run-specific COMBINED haplotype report into numeric 2-column format.

5.  Can join the species group info back to each individual in R after converting to 2-column format. The file “SebSppID\_3ChrCode.csv” will have all of the species for the individuals in the baseline (make sure to use file in the 11.10.2017 r project directory), and will just need to be modified to include the samples from your GTseq run of interest.

6.  Can then write this file to .csv in R. When open in excel, change all “NAs” in loci columns to “0” and remove the “.1” and “.2” after each locus in the header.

7.  Remove the 6 loci that were most prone to missing data and/or often have more than 3 haplotypes passing filter across individuals:

tag\_id\_1166 (often &gt;2 haplotypes) tag\_id\_934 (missing data and sometimes &gt;2 haplotypes) tag\_id\_2513 (appears to be some sort of repetitive element for some species outside of KGBC) tag\_id\_1871 (failed in many individuals) tag\_id\_1399 (failed in many individuals) tag\_id\_914 (failed in many individuals)

Then remove individuals missing genotype data at more than 15 loci. After converting to 2-column format and filtering there should be 1403 individuals in the baseline.

1.  This file will now include all baseline individuals as well as your (as yet unidentified samples) from your GTseq run. To run through gsi\_sim/rubias, will need to split this file into your “baseline” which is just the individuals used in creating the species ID vcf. The second file will be your “mixture” file, which just includes the samples from your GTseq run that you would like to assign using the baseline.

2.  Can then modify these files to match format required for whatever analyses you are interested in running (i.e. gsi\_sim, structure, etc.).
