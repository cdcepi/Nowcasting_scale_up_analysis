## A roadmap to account for reporting delays for public health situational awareness – a case study with COVID-19 and dengue in United States jurisdictions
This repository contains R scripts for all analyses included in the manuscript titled “A road map to account for reporting delays for public health situational awareness – a case study with COVID-19 and dengue in United States jurisdictions”. This work seeks to assess nowcast performance for varying epidemiologic surveillance systems to guide implementation across diverse systems. Our analytical goals are to 1) identify an optimized process for nowcasting when historical data are limited and reporting delays vary; 2) better understand and identify scenarios in which nowcasts fail; and 3) identify trade-offs in different nowcasting models. To meet these goals, this analysis used in 321 sequentially expanding datasets (n=280 for COVID-19 in the 6 states and n=41 for dengue in Puerto Rico) with unique ending dates.

### Repository structure
All included code can be run in R, with the required packages indicated at the top of each script. Analyses were performed in RStudio (RStudio 2023.12.0+369) using R (version 4.4.0) on a 64-bit Windows 10 desktop (128 GB RAM, Intel(R) Xeon(R) w5-3433 with 1.99 GHz processor).

All R scripts are in the **Code** folder and figures included in the manuscript are in the **Figures** folder. 

### Data sources
COVID-19 case data are available upon request at https://data.cdc.gov/Case-Surveillance/COVID-19-Case-Surveillance-Restricted-Access-Detai/mbd7-r32t

Dengue case data for 2010 are available from https://github.com/sarahhbellum/NobBS/blob/master/data/denguedat.RData and dengue case data from 2011-2014 were obtained from the Puerto Rico Department of Health.
