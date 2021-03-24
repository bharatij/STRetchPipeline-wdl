FROM ubuntu:16.04
MAINTAINER Bharati Jadhav <bharati.jadhav@mssm.edu>
LABEL \
  version="0.4.1" \
  description="STRetch image with modification of original estimateSTR.py script to perform outlier test on thousands of WGS and can be run on each chrosmosome separately"

# Update necessary packages
RUN apt-get update && apt-get install -qqy \
      git \
      default-jre-headless \
      zlib1g-dev \
      libncursesw5-dev \
      liblzma-dev \
      libbz2-dev \
      vim \
      ca-certificates \
      wget \
      make \
      bzip2 \
      unzip
     
#RUN install_packages default-jdk
RUN apt-get clean all \
	&& rm -rf /var/lib/apt/lists/* 
	
# based on https://github.com/Oshlack/STRetch/wiki/Installing-STRetch

# STRetch 

RUN wget "https://github.com/bharatij/STRetch/archive/master.zip" \
    && unzip master.zip \
    && rm master.zip

WORKDIR STRetch-master

#RUN ln -s /STRetch-master /STRetch

# prevent install.sh from downloading the 8Gb reference data zip - this should be downloaded outside the docker image
RUN mkdir -p reference-data && touch reference-data/hg19.PCRfreeWGS_143_STRetch_controls.tsv

RUN ./install.sh

RUN rm -rf reference-data/
WORKDIR ..