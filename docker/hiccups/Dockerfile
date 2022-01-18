FROM nvidia/cuda:8.0-devel-ubuntu16.04

LABEL maintainer "Paul Sud"
LABEL maintainer.email "encode-help@lists.stanford.edu"

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && \
    apt-get install -y \
        bzip2 \
        curl \
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
    git checkout 7b21fd620ee1f07266206caa2a7992d08d51ba8e && \
    chmod +x CPU/* CPU/common/* misc/* && \
    find -mindepth 1 -maxdepth 1  -type d -not -name "CPU" -not -name ".git" -not -name "misc" | xargs rm -rf

# Install Juicer tools
RUN curl \
        -L \
        https://github.com/aidenlab/Juicebox/releases/download/v2.13.06/juicer_tools_2.13.06.jar \
        -o /opt/juicer/CPU/common/juicer_tools.jar && \
    chmod 666 /opt/juicer/CPU/common/juicer_tools.jar && \
    ln -s juicer/CPU scripts

RUN curl \
        -LO \
        https://github.com/aidenlab/Juicebox/releases/download/v.2.14.00/feature_tools.jar

# For sorting, LC_ALL is C
ENV LC_ALL C
ENV PATH=/opt:/opt/scripts:/opt/scripts/common:$PATH
