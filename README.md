SRCommunityEpi_2019
====================

This repository contains code used in the following paper.

Quirine A. ten Bosch,Joseph M. Wagman, Fanny Castro-Llanos, Nicole L. Achee, John P. Grieco, T. Alex Perkins (2019) **Community-level impacts of spatial repellents for control of diseases vectored by Aedes aegypti mosquitoes**. *BioRXivs* 

All code contained within this repository is released under the [CRAPL v0.1 License](http://matt.might.net/articles/crapl/). Data provided in this repository are sufficient to rerun all analyses. This code base is part of an associated Open Science Framework project [https://osf.io/j9cks/](https://osf.io/j9cks/).

The scripts for simulations and model fitting were performed on the University of Notre Dame's Center for Research Computing cluster (experimental huts studies) [http://crc.nd.edu](http://crc.nd.edu) and on a desktop computer (Mac OSX). Processing of outputs was also done on a desktop computer (Mac OSX) 

====================

### Set up of code base: 

* `Code` contains code to analyze the data, perform calculations, and generate output and figures 
* `Data` contains experimental data described in manuscripts. Read by scripts from the `Code` folder
* `Input` contains RData files with workspaces used for final 'FOI' portion of the analyses. 

### Code folder

The scripts in this folder are used to analyze experimental data and use estimated parameters in FOI-framework. 

Bloodfeeding:
* `Bloodfeeding.R`, which read and analyzes data from `./Data/bloodfeeding.csv`

Longevity: 
* `LongevityResidual.R`, which calls `Longevity_PrepareData.R` to read in and process data and then performs conditional parameteric survival analyses 
* `LongevityInstantaneous.R`, which calls `Longevity_PrepareData.R` to read in and process data and then estimates instantaneous death rates by dosage 
* `Longevity_PrepareData.R`, reads in data from `./Data/Longevity_XXFAR.csv`, converts these to a time-to-event format, and combines different dosages into one data.frame
* `LongevityFunctions.R`, is called by other Longevity scripts and contains functions for data processing and analysis 

Movement:
* `Movement.R`, reads in RData file from `./Input/Movement.RData` and generates summary statistics and figures (SFig2)
Other code and data behind these numbers are deposited elsewhere [https://osf.io/5hcpf/](https://osf.io/5hcpf/)

Force of Infection: 
* `RelFOI.R` , loads data from movement, longevity, and bloodfeeding experiments (stored in `./Input`) and performs all FOI-relatved calculations and creates FOI-figures (3,4,S4)

### Data folder
* `SR_Bloodfeeding.csv` , contains data from bloodfeeding experiments: with the numbers fully, partially, and non-fed by time and dosage
* `Longevity_XXFAR.csv` , contains data from longevity experiments: number of mosquitoes surviving over time for up to 6 replicates of control and treatment arms
Data from movement experiments can be found on Open Science Framework [https://osf.io/5hcpf/](https://osf.io/5hcpf/)

### Input folder 
* `BloodFeeding.RData`, data generated by `Bloodfeeding.R`
* `Lethality.Instantaneous.RData`, data generated by `LongevityResidual.R`
* `Lethality.Residual.Exp.RData`, data generated by `LongevityInstantaneous.R`
* `Movement.RData`, data generated as part of manuscript on experimental hut experiments
Quirine A. ten Bosch,Fanny Castro-Llanos, Hortance Manda, Amy C. Morrison, John P. Grieco, Nicole L. Achee, T. Alex Perkins (2018) **Model-based analysis of experimental hut data elucidates multifaceted effects of a volatile chemical on Aedes aegypti mosquitoes**. *Parasites and Vectors* [ https://doi.org/10.1186/s13071-018-2919-0]( https://doi.org/10.1186/s13071-018-2919-0) 
with code available on [https://github.com/quirine/ExperimentalHuts](https://github.com/quirine/ExperimentalHuts) 
and data on [https://osf.io/5hcpf/](https://osf.io/5hcpf/)

Any questions about this code base can be directed to qtenbosc@alumni.nd.edu
