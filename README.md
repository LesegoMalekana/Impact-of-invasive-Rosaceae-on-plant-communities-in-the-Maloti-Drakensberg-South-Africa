Data, figures, and R scripts for the manuscript examining the impact of invasive Rosaceae species on plant communities and microclimate in the Maloti-Drakensberg, South Africa.

Project description

This repository contains the datasets and R scripts used to analyse the effects of invasive Rosaceae species on plant diversity, species composition, vegetation cover, and microclimate in montane grasslands of the Maloti-Drakensberg, South Africa.

Data
**01_Diversity_Impact_Analysis_and_NMDS**

Dataset used to analyse diversity metrics and community composition (NMDS).

Variables

Elevation: Elevation category (three levels)
Site: Study site identifier
Date: Sampling date
Treatment: Invasion status (Invaded or Control)
Plot_ID: Plot identifier
Natives_sp: List of native species recorded in each plot
Natives_nr: Total abundance (number of individuals) of native species
Alien_sp: List of non-native species recorded in each plot
Alien_nr: Total abundance (number of individuals) of non-native species

**02_Cleaned_Light_Temperature_Standardised_Publish**

Dataset containing daily light and temperature measurements from invaded and control plots.

Variables

Site: Study site identifier
Treatment: Invasion status (Invaded or Control)
Datetime: Date and time of measurement
Temperature_c: Air temperature in degrees Celsius
Light_lux: Light intensity in lux

**03_Cover_Change_Perc**

Dataset containing percentage cover of vegetation groups across treatments.

Variables

Site: Study site identifier
Date: Sampling date
Treatment: Invasion status (Invaded or Control)
Native: Native vegetation cover
Alien: Non-native vegetation cover
Focal: Cover of focal Rosaceae species
Grasses: Graminoid cover
Bare_ground: Bare ground cover
R scripts

**01_Diversity_Impact_Analysis_and_NMDS.R**

Script used to analyse diversity metrics, beta diversity, and community composition (NMDS).

**02_Light_Temperature_Analysis.R**

Script used to analyse differences in light intensity and air temperature between invaded and control plots.

**03_Cover_Change_Veg_Bare.R**

Script used to visualise differences in vegetation cover and bare ground between invaded and control plots.

**Figures**

This folder contains the figures used in the manuscript. Most figures are generated automatically when the scripts are run.

**Tables**

Statistical tables are printed to the R console during analysis. The scripts also include options to export results to Excel files.

**Notes on data processing**
All datasets were cleaned prior to upload.
Missing values are represented as NA.

**Reproducibility**

All analyses were conducted in R. The scripts provided in this repository reproduce the main results presented in the manuscript. Users can download the data and scripts and run the analyses locally. File paths in the scripts may need to be adjusted to match the user’s local directory structure.
