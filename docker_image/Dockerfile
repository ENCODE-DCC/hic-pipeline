############################################################
# Dockerfile to build Juicer container image
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu:16.04 as main

LABEL maintainer "Paul Sud"
LABEL maintainer.email "encode-help@lists.stanford.edu"

# Update the repository sources list
# Install base packages: java, git, wget
RUN apt-get update && apt-get install -y \
    openjdk-8-jdk \
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
RUN git clone --branch encode https://github.com/theaidenlab/juicer.git && \
    cd juicer && \
    git checkout 44752230566dffae0aabddb35116f5d9b72aface && \
    chmod +x CPU/* CPU/common/*

# Install Juicer tools
ADD http://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.11.09_jcuda.0.8.jar /opt/juicer/CPU/common
RUN ln -s /opt/juicer/CPU/common/juicer_tools_1.11.09_jcuda.0.8.jar /opt/juicer/CPU/common/juicer_tools.jar
RUN chmod 666 /opt/juicer/CPU/common/juicer_tools_1.11.09_jcuda.0.8.jar
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

# Use multistage build for fragment
FROM rust:1.42.0-slim-stretch@sha256:1201e100253c82c4accd19c649d187ab741b89954397dcdddb792008f898f93a as builder
WORKDIR /hic-pipeline
COPY . .
RUN cargo build --release

FROM main

RUN mkdir -p hic-pipeline/hic_pipeline
COPY hic_pipeline hic-pipeline/hic_pipeline/
# Add compiled binary from other stage of build
COPY --from=builder /wgbs-pipeline/target/release/fragment hic-pipeline/hic_pipeline/
ENV PATH="/opt/hic-pipeline/hic_pipeline:${PATH}"
