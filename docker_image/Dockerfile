############################################################
# Dockerfile to build Juicer container image
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu:16.04

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
    python3 \
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
ADD http://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.11.04_jcuda.0.8.jar /opt/juicer/CPU/common
RUN ln -s /opt/juicer/CPU/common/juicer_tools_1.11.04_jcuda.0.8.jar /opt/juicer/CPU/common/juicer_tools.jar
RUN chmod 666 /opt/juicer/CPU/common/juicer_tools_1.11.04_jcuda.0.8.jar
RUN ln -s juicer/CPU scripts

# Install Straw for removal of hic header in create_hic test
RUN git clone --branch master https://github.com/aidenlab/straw.git

# Install pairix
ADD https://github.com/4dn-dcic/pairix/archive/0.3.6.zip .
RUN mv 0.3.6.zip pairix-0.3.6.zip && unzip pairix-0.3.6.zip 
RUN cd pairix-0.3.6/ && make

RUN rm pairix-0.3.6.zip && rm samtools-1.6.tar && rm v0.7.17.zip

# For sorting, LC_ALL is C
ENV LC_ALL C
ENV PATH=/opt:/opt/scripts:/opt/scripts/common:$PATH

# Largely copied from https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docker_image/Dockerfile
# Get scripts from hic-pipeline repo
# This COPY assumes the build context is the root of the hic-pipeline repo
# and it gets whatever is checked out plus local modifications
# so the buildling command should:
# cd [GIT_REPO_DIR] && docker build -f docker_image/Dockerfile .
RUN mkdir -p hic-pipeline/src
COPY src hic-pipeline/src/
ENV PATH="/opt/hic-pipeline/src:${PATH}"

#CMD ["-h"]


