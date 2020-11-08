# DOCKERFILE for GEX
# github.com/davidsbutcher/gex
# docker build --tag davidsbutcher/gex:releaseX --no-cache .

# get r-base image
FROM rocker/r-base:4.0.3

#### install R packages from github

# RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"

RUN R -e "devtools::install_github('davidsbutcher/GEX')"

#### install R packages from Bioconductor

# RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"

#### install R packages from CRAN 

RUN R -e "install.packages('magrittr', repos='http://cran.rstudio.com/')"

#### copy the app to the image
COPY ui.R /var/gex/
COPY server.R /var/gex/
COPY about.md /var/gex/
COPY /www/ /var/gex/www/

#### select port
EXPOSE 3838

#### run app
CMD R -e 'shiny::runApp("/var/gex", port = 3838, host = "0.0.0.0")'