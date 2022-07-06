# Gossip in Hungarian firms

This repository contains the replication package for the article: 
Estévez, J. L., & Takács, K. (2022). Brokering or sitting between two chairs? A group perspective on workplace gossip. *Frontiers in Psychology, 13*: 815383. https://doi.org/10.3389/fpsyg.2022.815383.

**Software requirements**

R (code was last run with R version 4.1.2 in RStudio 2021.09.0)
- ggplot2 (3.3.5)
- lpSolve (5.6.15)
- irr (0.84.1)
- ape (5.6-1)
- lattice (0.20-44)
- viridis (0.6.2)
- sna (2.6)
- igraph (1.2.11)
- ggpubr (0.4.0)
- lme4 (1.1-27.1)
- effectsize (0.6.0.1)
- insight (0.17.1)

**File list**

- 1_Data_tidying.R
- 2_Composite_networks.R
- 3_Descriptive_analysis.R
- 4_ML_analysis.R
- data.RData

**Instructions**

Open RStudio and create a new project, preferably in a new directory (the name of this project is up to you). This creates a new folder with a single .Rproj file named after the name you chose. In this folder, drag all the files mentioned in the file list. After this, both the data and code must be available within RStudio (see tab “files”) and can be accessed by simply clicking on them.

The code needs to be run in order. The first file (1_Data_tidying.R) will prepare the data for descriptive analyses. After running this file, a new .RData file (tidieddata.RData) will appear in the folder. This file is the input for 2_Composite_networks.R. After running 2_Composite_networks.R, the file tidieddata2.RData is created, which is the input for 3_Descriptive_analysis.R. Finally, running 3_Descriptive analysis.R creates modellingdata.R, the input for the last file: 4_ML_analysis.

For the creation of composite networks, run 2_Composite_networks.R

For descriptive analyses, run  3_Descriptive_analysis.R

For the main analyses, run 4_ML_analyis.R
