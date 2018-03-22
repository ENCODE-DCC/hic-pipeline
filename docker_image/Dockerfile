############################################################
# Dockerfile to build Juicer container image
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu

# File Author / Maintainer
MAINTAINER Neva C. Durand 

# Update the repository sources list
# Install base packages: java, git, wget
RUN apt-get update && apt-get install -y \
    default-jdk \
    gawk \
    gcc \
    git \
    libz-dev \
    locales \
    make \
    unzip \
    bzip2 \
    libbz2-dev \	
&& rm -rf /var/lib/apt/lists/*

# GAWK has the 'and' function, needed for chimeric_blacklist
RUN echo 'alias awk=gawk' >> ~/.bashrc

# Need to be sure we have this for stats
RUN locale-gen en_US.UTF-8

WORKDIR /opt/

# Install BWA
ADD https://github.com/lh3/bwa/archive/v0.7.17.zip .
RUN unzip v0.7.17.zip 
RUN cd bwa-0.7.17/ && make
RUN ln -s bwa-0.7.17/bwa bwa

# Install Samtools
ADD https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 .
RUN bunzip2 samtools-1.6.tar.bz2 
RUN tar xf samtools-1.6.tar 
RUN cd samtools-1.6 && ./configure --without-curses --disable-bz2 --disable-lzma && make && make install

# Install Juicer
RUN git clone --branch encode https://github.com/theaidenlab/juicer.git
RUN cd juicer && chmod +x CPU/* CPU/common/* 

# Install Juicer tools
ADD http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.7.6_jcuda.0.8.jar /opt/juicer/CPU/common
RUN ln -s /opt/juicer/CPU/common/juicer_tools.1.7.6_jcuda.0.8.jar /opt/juicer/CPU/common/juicer_tools.jar
RUN ln -s juicer/CPU scripts

# Version number contained in image
ADD VERSION .

# For sorting, LC_ALL is C
ENV LC_ALL C
ENV PATH=/opt:/opt/scripts:/opt/scripts/common:$PATH

ENTRYPOINT ["juicer.sh", "-D", "/opt"]
CMD ["-h"]


