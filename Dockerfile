FROM ubuntu:16.04

MAINTAINER Kristian Peters (kpeters@ipb-halle.de)

# Add cran R backport
ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get -y update && apt-get -y dist-upgrade && apt-get -y install apt-transport-https perl
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Generate locales
ENV LC_ALL="en_US.UTF-8"
ENV LC_CTYPE="en_US.UTF-8"
RUN locale-gen $LC_ALL
RUN dpkg-reconfigure locales

# Install packages
RUN apt-get -y update && apt-get -y dist-upgrade && apt-get -y install apt-transport-https make gcc gfortran g++ libblas-dev liblapack-dev libxml++2.6-dev libexpat1-dev libxml2-dev libnetcdf-dev libssl-dev r-base r-base-dev maven texlive-latex-base git openjdk-8-jdk-headless texlive-fonts-recommended openjdk-8-jre-headless pkg-config parallel wget curl git unzip zip

# Install R packages
RUN R -e 'install.packages(c("irlba","igraph","ggplot2","digest","lattice","XML","Rcpp","reshape2","plyr","stringi","stringr","tools","intervals","devtools","RColorBrewer","plyr","RANN","knitr","ncdf4","microbenchmark","RUnit","parallel","foreach","doMC","curl","jsonlite"), repos="https://cloud.r-project.org/")'

# Install  Bioconductor 
RUN R -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("multtest","MSnbase","mzR","MassSpecWavelet","S4Vectors","BiocStyle","faahKO","msdata","xcms","CAMERA"), ask=FALSE)'

# Install MetFrag
RUN wget http://central.maven.org/maven2/net/sf/jni-inchi/jni-inchi/0.8/jni-inchi-0.8.jar && mkdir -p /root/.jnati/repo/ && jar xf jni-inchi-0.8.jar && mv META-INF/jniinchi /root/.jnati/repo/
RUN wget -O /usr/local/bin/MetFragCLI.jar http://msbi.ipb-halle.de/~cruttkie/92f73acb731145c73ffa3dfb8fd59581bee0d844963889338c3ec173874b5a5f/MetFrag-2.4.5.jar
RUN wget -O /usr/local/bin/metfrag.jar https://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.5-CL.jar
RUN wget -O /usr/local/bin/metfrag https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/metfrag/metfrag.sh

# Install MetFrag tools
RUN git clone https://github.com/c-ruttkies/ConvertMetFragCSV.git && cd ConvertMetFragCSV && mvn clean install -am && mv target/ConvertMetFragCSV-*-jar-with-dependencies.jar /usr/local/bin/ConvertMetFragCSV.jar 
RUN wget -O /usr/local/bin/MetFrag-Tools.jar https://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.5-Tools.jar

# Install MetFrag lib
RUN wget -O /usr/local/bin/MetFrag-Lib.jar https://msbi.ipb-halle.de/~cruttkie/metfrag/MetFragLib-2.4.5.jar

# Install metfrag-cli-batch
RUN wget -O /usr/local/bin/run_metfrag.sh https://raw.githubusercontent.com/phnmnl/container-metfrag-cli-batch/develop/run_metfrag.sh

# Cleanup
RUN apt-get -y --purge --auto-remove remove make gcc gfortran g++
RUN apt-get -y clean && apt-get -y autoremove && rm -rf /var/lib/{cache,log}/ /tmp/* /var/tmp/*

# Add scripts folder to container
ADD scripts/* /usr/local/bin/
RUN chmod +x /usr/local/bin/*

