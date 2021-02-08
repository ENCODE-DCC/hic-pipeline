FROM nvidia/cuda:8.0-devel-ubuntu16.04

LABEL maintainer "Paul Sud"
LABEL maintainer.email "encode-help@lists.stanford.edu"

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && \
    apt-get install -y \
        bzip2 \
        gawk \
        gcc \
        git \
        libbz2-dev \
        libz-dev \
        locales \
        make \
        openjdk-11-jdk \
        unzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/

# Install Juicer
RUN git clone --branch encode https://github.com/theaidenlab/juicer.git && \
    cd juicer && \
    git checkout 7d482a7dfb261ac6193e653144725958d483ec40 && \
    chmod +x CPU/* CPU/common/* misc/* && \
    find -mindepth 1 -maxdepth 1  -type d -not -name "CPU" -not -name ".git" -not -name "misc" | xargs rm -rf

# Install Juicer tools
ADD juicer_tools.2.06.03.jar /opt/juicer/CPU/common/juicer_tools.jar
RUN chmod 666 /opt/juicer/CPU/common/juicer_tools.jar
RUN ln -s juicer/CPU scripts

# Version number contained in image
# ADD VERSION .

# For sorting, LC_ALL is C
ENV LC_ALL C
ENV PATH=/opt:/opt/scripts:/opt/scripts/common:$PATH
