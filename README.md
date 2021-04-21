# LOCeating_AssociationRulesMining
The project used data compiled across multiple studies to identify characteristics associated with loss of control (LOC)-eating in children. A data-driven approach called association rules mining (ARM) was used to generate clusters of characteristics that are associated with LOC-eating. This study served to generate new hypotheses about LOC-eating in children.

    Copyright (C) 2020 Alaina L Pearce

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


## Primary Publication
These data have now been published in Appetite:
Pearce, A. L., Brick, T. R., Masterson, T., Adise, S., Fearnbach, S. N., Stein, W., English, L., Tanofsky-Kraff, M., & Keller, K. L. (2021). Using association rules mining to characterize loss of control eating in childhood. Appetite, 163, 105236. https://doi.org/10.1016/j.appet.2021.105236
*While there is an embargo on this paper, a pre-proof version will be available on my personal website: https://sites.psu.edu/alainapearce/research/publications/

A perminant link for this publications data and procedures is register through Open Science Framework (10.17605/OSF.IO/BJ8ZA) with web link: https://osf.io/bj8za/. Additionally, the version of code used for this publication and submitted to OSF is stored in the file 2020_LOC_ARM_OSF directory in this repository. While I will not edit anything in that directory, I cannot garentee other docuements in this repository will not be updated for future analyses. 

## Disclaimers
These scripts are not guaranteed to work for new data or under new directory configurations, however, they should work if  downloaded with the data and no changes are made to directories. To use these scripts R/Rstudio and all required libraries. I complete this project prior to learning abour renv so exact replication of results is not gaurenteed as both R and package versions may differ.

## Directory Structure
-2020_LOC_ARM_OSF: contains exact version of code and data used for the Pearce et al. (2020) publication
-Data: contains raw data
-ResultsOutput: .csv files for the resulting rules
-LOC_Figs: figures generated by the .Rmd file
-Primary Directory: contains all scripts
  1) LOCeating_ARM.Rproj - the Project file used to analyze data
  2) LOC_ARules_setup_CascadeApproach.R - processes and cleans data; also creates a dataset with no duplicates (takes first visit only if there were multiple)
  3) LOC_Analyses_CascadeApproach.R - primary script; contains all the analysis code
  4) LOC_ARules_CascadeApproach.Rmd - markdown file that will source the analysis script and produce a .pdf summary of results
  5) functions.R - a misc. set of custom functions; not all are used in this project
  6) ARM_ORconf.R - a script to calcuate the confidence intervals of the odds ratios
