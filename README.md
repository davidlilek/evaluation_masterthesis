# Data analysis and evaluation master thesis

This repo contains all data evaluation steps which were performed in the course of the master thesis "Development and validation of a workflow for the evaluation of bottom-up
proteomics data using MaxQuant and R" in 2022. The master thesis was submitted by David Lilek at the Austrian Biotech University of Applied Sciences.

Most of the customized scripts were created using R for post-processing of the MaxQuant results. The post-processing included filterung, statistical analysis and vizualization.
R version R version 4.2.0 was used ([session_info_R.txt](https://github.com/davidlilek/evaluation_masterthesis/files/9149696/session_info_R.txt)).

For performance test in MaxQuant different parameters were varied, as shown below. Derived from this result in the script names. For example:

* MBR_pooled_oldversion: MBR was used; pooling was done in MaxQuant; old version of MaxQuant was used
* MBR: MBR was used; pooling was done manually
* noMBR_pooled_secondpeptides: MBR was disabled; pooling was done in MaxQuant; option second peptides was enabled
* ...

<img width="517" alt="performancetest_maxquant_pic" src="https://user-images.githubusercontent.com/60740660/179982995-e85ee211-f597-4a03-9692-1d3d102398fa.png">

The different evaluation methods were analyzed using different protein types as shown below. 

* Razor: non-unique peptides; these peptides are assigned to the protein group with the largest number of total identified peptides |
* Unique: peptides which are unique for this protein group | 
* LFQ: label free quantification (quantified proteins). 

The option one peptide or two peptides refers to the minimum number of peptides required for protein identification.

<img width="374" alt="proteintypes_github" src="https://user-images.githubusercontent.com/60740660/180036892-1f627191-42c4-4117-9e3d-fd464d4f8538.png">

## QC

Contains one Rmd for each evaluation method as described above resulting in 8 Rmd files. `20220406_TR_QC` refers to the LC-MS/MS sequence name.
The results are summarized in 4 Rmd files.

* `20220406_TR_QC_summary_rerun_onepeptide.Rmd`
* `20220406_TR_QC_summary_rerun_seconpeptides_onepeptide.Rmd`
* `20220406_TR_QC_summary_rerun_twopeptide.Rmd`
* `20220406_TR_QC_summary_seconpeptides_rerun_twopeptide.Rmd`

In this files the results are sumarized and plots created and arranged for master thesis in file `afterprocessing4masterthesis.R`.

LFQ intensities were compared and analyzed using `compare_lfqintensities.R`.

## UPS1

Contains one Rmd for each evaluation method as described above resulting in 8 Rmd files. `20220406_TR_QC` refers to the LC-MS/MS sequence name.
The results are summarized in 4 Rmd files.

* `UPS1_summary_rerun_secondpeptides_onepeptide.Rmd`
* `UPS1_summary_rerun_secondpeptides_onepeptide.Rmd`
* `UPS1_summary_rerun_twopeptide.Rmd`
* `UPS1_summary_rerun_secondpeptides_twopeptide.Rmd`

In this files the results are sumarized and plots created.

LFQ intensities were compared and analyzed using `compare_lfqintensities.R`.

All results are summarized and visualized in `afterprocessing4masterthesis.R`. 

* single plots for each concentration
* analysis and comparison using one and two peptides for post-processing
* comparison number of human proteins
* analysis and comaprison of different evaluation methods
* comparison of MBR and no MBR

## Extracts

## scripts-4masterthesis




