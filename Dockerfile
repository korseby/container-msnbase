FROM container-registry.phenomenal-h2020.eu/phnmnl/camera:latest

MAINTAINER PhenoMeNal-H2020 Project (phenomenal-h2020-users@googlegroups.com)

LABEL software=MSnbase
LABEL software.version=2.0.2
LABEL version=0.1
LABEL Description="MSnbase: Basic plotting, data manipulation and processing of MS-based Proteomics data."

# Install packages for compilation
RUN apt-get -y update
RUN apt-get -y --no-install-recommends install make gcc gfortran g++ libblas-dev liblapack-dev libxml++2.6-dev libexpat1-dev

# Install dependencies
RUN R -e 'install.packages(c("ggplot2","digest","lattice","XML","Rcpp","reshape2","plyr","stringr","intervals"), repos="https://mirrors.ebi.ac.uk/CRAN/")'

# Install CAMERA
RUN R -e 'source("https://bioconductor.org/biocLite.R"); biocLite("MSnbase")'

# De-install not needed packages
RUN apt-get -y --purge --auto-remove remove make gcc gfortran g++

# Clean-up
RUN apt-get -y clean && apt-get -y autoremove && rm -rf /var/lib/{cache,log}/ /tmp/* /var/tmp/*

# Add scripts folder to container
#ADD show_chromatogram.r /usr/local/bin/show_chromatogram.r
#RUN chmod +x /usr/local/bin/show_chromatogram.r

# Define Entry point script
#ENTRYPOINT [ "Rscript" ]
#CMD [ "/usr/local/bin/show_chromatogram.r" ]
