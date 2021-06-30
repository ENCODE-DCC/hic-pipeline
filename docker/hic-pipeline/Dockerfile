FROM ubuntu:16.04@sha256:e9938f45e51d9ff46e2b05a62e0546d0f07489b7f22fbc5288defe760599e38a as main

LABEL maintainer "Paul Sud"
LABEL maintainer.email "encode-help@lists.stanford.edu"

# Package versions
ENV BWA_VERSION=0.7.17
ENV PAIRIX_VERSION=0.3.6
ENV SAMTOOLS_VERSION=1.6

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && \
    apt-get install -y \
        bzip2 \
        curl \
        g++ \
        gawk \
        gcc \
        git \
        libbz2-dev \
        libcurl4-openssl-dev \
        libz-dev \
        locales \
        make \
        openjdk-11-jdk \
        python3 \
        unzip \
    && rm -rf /var/lib/apt/lists/*

# GAWK has the 'and' function, needed for chimeric_blacklist
RUN echo 'alias awk=gawk' >> ~/.bashrc

# Need to be sure we have this for stats
RUN locale-gen en_US.UTF-8

WORKDIR /opt/

# Install BWA
RUN curl -OL "https://github.com/lh3/bwa/archive/v${BWA_VERSION}.zip" && \
    unzip "v${BWA_VERSION}.zip" && \
    cd "bwa-${BWA_VERSION}/" && \
    make && \
    cp bwa /usr/local/bin && \
    cd .. && \
    rm -rf "bwa-v${BWA_VERSION}"

# Install Samtools
RUN curl -OL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    bunzip2 "samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    tar xf "samtools-${SAMTOOLS_VERSION}.tar" && \
    cd "samtools-${SAMTOOLS_VERSION}" && \
    ./configure --without-curses --disable-bz2 --disable-lzma && \
    make && \
    make install && \
    cd .. && \
    rm -rf "samtools-${SAMTOOLS_VERSION}" "samtools-${SAMTOOLS_VERSION}.tar"

RUN curl -OL "https://github.com/4dn-dcic/pairix/archive/${PAIRIX_VERSION}.zip" && \
    unzip "${PAIRIX_VERSION}.zip" && \
    cd "pairix-${PAIRIX_VERSION}/" && \
    make && \
    cp bin/pairix bin/bgzip /usr/local/bin && \
    cp util/merged_nodup2pairs.pl /opt && \
    cd .. && \
    rm -rf "${PAIRIX_VERSION}.zip" "pairix-${PAIRIX_VERSION}"

# Eigenvector
RUN git clone https://github.com/moshe-olshansky/EigenVector.git && \
    cd EigenVector && \
    git checkout 68c399f4ac7d1b953724ab41ae0c4b32368392b1 && \
    git clone https://github.com/aidenlab/straw.git && \
    cd straw && \
    git checkout 7d5ec096e116d6680455bcc81eee8e42c2adf6ba && \
    cp C++/straw.cpp C++/straw.h ../C++/FlipSign && \
    cd ../C++/FlipSign && \
    g++ \
        -O \
        --std=c++11 \
        -o newGW_Intra_Flip \
        GWevIntra_new.cpp \
        theBestEigen.c \
        thdMul.c \
        hgFlipSign.c \
        straw.cpp \
        -I. \
        -lz \
        -lcurl \
        -lpthread && \
    chmod +x newGW_Intra_Flip  && \
    mv newGW_Intra_Flip  /usr/local/bin && \
    rm -rf newGW_Intra_Flip

RUN curl \
    -OL \
    https://github.com/ENCODE-DCC/kentUtils_bin_v381/raw/d0ba32f3aee6da02c0f75d8f998c5f1d921b047c/bin/wigToBigWig && \
    chmod +x wigToBigWig && \
    mv wigToBigWig /usr/local/bin

# Install Juicer
RUN git clone --branch encode https://github.com/theaidenlab/juicer.git && \
    cd juicer && \
    git checkout 3cbd828fff8d2ae9ac675d517fa8139038167123 && \
    chmod +x CPU/* CPU/common/* misc/* && \
    find -mindepth 1 -maxdepth 1  -type d -not -name "CPU" -not -name ".git" -not -name "misc" | xargs rm -rf

# Install Juicer tools
RUN curl \
        -L \
        https://github.com/aidenlab/Juicebox/releases/download/v2.13.06/juicer_tools_2.13.06.jar \
        -o /opt/juicer/CPU/common/juicer_tools.jar && \
    chmod 666 /opt/juicer/CPU/common/juicer_tools.jar && \
    ln -s juicer/CPU scripts

# For sorting, LC_ALL is C
ENV LC_ALL C
ENV PATH=/opt:/opt/scripts:/opt/scripts/common:/opt/juicer/misc:$PATH

RUN mkdir -p hic-pipeline/hic_pipeline
COPY hic_pipeline hic-pipeline/hic_pipeline/
ENV PATH="/opt/hic-pipeline/hic_pipeline:${PATH}"
