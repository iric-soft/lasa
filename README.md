This is the web interface to LASA: Leucegene AML Surfaceome Atlas.

It's a R Shiny app which launches severals analyses on Leucegene surfaceome cohort.

## Cite

 - Paper: [Immunotherapeutic targeting of surfaceome heterogeneity in AML]()
 - DOI of the code using, Zenodo: [![DOI](https://zenodo.org/badge/739109138.svg)](https://zenodo.org/doi/10.5281/zenodo.10460001)

## Requirement

 - R version 4.0.0 (DEP, used to analyse surfaceome data, doesn't work on R > 4.0.0)

## Data

All data used in LASA can be found in the **Data** folder, except for single cell data which are too large. To get single cell data check the download section of [LASA](https://lasa.leucegene.ca) and copy/extract the downloaded file **as** data/singleCell/sc_samples.loom

**Data folder will be added when the paper will be accepted.**

For more details or to generate data by your self, see the method section of [Immunotherapeutic targeting of surfaceome heterogeneity in AML]().

## Install

```
# Create empty folders, useful to launch LASA using a shiny server
mkdir -p shiny/bookmark
mkdir -p shiny/log

# Edit R_LIBS path in .Renviron
nano shiny/.Renviron

# Copy .Renviron in severals place (to simplify your life)
ln -s ./shiny/.Renviron ./.Renviron
cd ./shiny/www/lasa
ln -s ../../.Renviron ./.Renviron
cd -

# Create empty folder and install R dependencies
mkdir -p shiny/R_lib
cd shiny/
Rscript install_packages.r
cd -
```


## Launch (locally) using rstudio
```
# In main directory launch rstudio
rstudio &
```

Open `shiny/www/lasa/server.R` file and click on "run App" button.

## Launch using a shiny server

```
cd shiny

# Edit/Create a config file in etc folder
cp etc/shiny-server.conf etc/shiny-server.lasa.conf
nano etc/shiny-server.lasa.conf

# Run shiny server
nohup /path/to/shiny-server etc/shiny-server.lasa.conf < /dev/null &>> ./log/shiny_lasa.log &

```
