# Impact-of-invasive-Rosaceae-on-plant-communities-in-the-Maloti-Drakensberg-South-Africa
Data file and R code for the paper


Project description
This repository contains data and R scripts used to analyse the effects of invasive Rosaceae species on vegetation composition and microclimate in montane grasslands.

Data files

1. Cleaned_Light_Temperature_Standardised_Publish.csv
Cleaned dataset containing light intensity and temperature measurements at the plot level.

Variables:

Site: study site identifier
Treatment: invasion status (e.g. invaded, uninvaded)
Datetime: date and time of measurement
Temperature_c: air temperature in degrees Celsius
Light_lux: light intensity in lux

2. Community_Composition_Publish.csv
Dataset containing species composition and abundance per plot.

Variables:

Elevation: elevation category (three levels)
Site: study site identifier
Date: sampling date
Treatment: invasion status (e.g. invaded, uninvaded)
Plot_ID: plot identifier
Natives_sp: list of native species recorded in the plot
Natives_nr: total abundance (number of individuals) of native species
Alien_sp: list of non-native species recorded in the plot
Alien_nr: total abundance (number of individuals) of non-native species

Note: Species lists are stored at the plot level. Abundance values represent counts of individuals per plot.

R scripts

1. Impact_of_Rosaceae_Trees_Publish.R
Script used to analyse the effects of Rosaceae invasion on vegetation composition and abundance.

2. Light_and_Temperature_Analyses_Publish.R
Script used to analyse differences in light intensity and temperature between invaded and uninvaded plots.

Notes on data processing

All datasets have been cleaned prior to upload.
Missing values are represented as NA.

Reproducibility

All analyses were conducted in R. Scripts are provided to reproduce the main results presented in the manuscript. Users can download the data and scripts and run the analyses locally. File paths in the scripts may need to be adjusted depending on the user’s local directory structure.
