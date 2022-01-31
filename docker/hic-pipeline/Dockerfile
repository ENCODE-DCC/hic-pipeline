FROM ubuntu:16.04@sha256:e9938f45e51d9ff46e2b05a62e0546d0f07489b7f22fbc5288defe760599e38a as main

LABEL maintainer "Paul Sud"
LABEL maintainer.email "encode-help@lists.stanford.edu"

# Package versions
ENV BWA_VERSION=0.7.17
ENV PAIRIX_VERSION=0.3.6

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && \
    apt-get install -y \
        bc \
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
        parallel \
        python3 \
        unzip \
    && rm -rf /var/lib/apt/lists/*

# GAWK has the 'and' function, needed for chimeric_blacklist
RUN echo 'alias awk=gawk' >> ~/.bashrc

RUN ln -s /usr/bin/python3 /usr/bin/python

# Need to be sure we have this for stats
RUN locale-gen en_US.UTF-8

# Fix warning for diploidify
# https://www.educative.io/edpresso/error-mesg-ttyname-failed-inappropriate-ioctl-for-device
# Need to escape the &s
# https://unix.stackexchange.com/questions/32907/what-characters-do-i-need-to-escape-when-using-sed-in-a-sh-script
RUN sed -i 's/mesg n || true/tty -s \&\& mesg n/' /root/.profile

WORKDIR /opt/

# Install BWA
RUN curl -OL "https://github.com/lh3/bwa/archive/v${BWA_VERSION}.zip" && \
    unzip "v${BWA_VERSION}.zip" && \
    cd "bwa-${BWA_VERSION}/" && \
    make && \
    cp bwa /usr/local/bin && \
    cd .. && \
    rm -rf "bwa-v${BWA_VERSION}"

RUN export SAMTOOLS_VERSION=1.14 && \
    curl -OL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
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
    cp util/juicer_shortform2pairs.pl /opt && \
    cd .. && \
    rm -rf "${PAIRIX_VERSION}.zip" "pairix-${PAIRIX_VERSION}"

# Eigenvector
RUN git clone https://github.com/aidenlab/EigenVector.git && \
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
    mv newGW_Intra_Flip /usr/local/bin && \
    rm -rf newGW_Intra_Flip

RUN git clone https://github.com/ENCODE-DCC/kentUtils_bin_v381.git && \
    cd kentUtils_bin_v381/bin && \
    chmod +x wigToBigWig bedGraphToBigWig && \
    mv wigToBigWig bedGraphToBigWig /usr/local/bin && \
    cd ../../ && \
    rm -rf kentUtils_bin_v381

# Install Juicer
RUN git clone --branch encode https://github.com/theaidenlab/juicer.git && \
    cd juicer && \
    git checkout 50d557f1d4725a475071fce5975839602bd311e5 && \
    chmod +x CPU/* CPU/common/* misc/* && \
    find -mindepth 1 -maxdepth 1  -type d -not -name "CPU" -not -name ".git" -not -name "misc" | xargs rm -rf

# Install Juicer tools
RUN curl \
        -L \
        https://github.com/aidenlab/Juicebox/releases/download/v2.13.06/juicer_tools_2.13.06.jar \
        -o /opt/juicer/CPU/common/juicer_tools.jar && \
    chmod 666 /opt/juicer/CPU/common/juicer_tools.jar && \
    ln -s juicer/CPU scripts && \
    ln -s /opt/juicer/CPU/common/juicer_tools /opt/juicer/CPU/juicer_tools && \
    ln -s /opt/juicer/CPU/common/juicer_tools.jar /opt/juicer/CPU/juicer_tools.jar

RUN curl \
        -LO \
        https://github.com/aidenlab/Juicebox/releases/download/v.2.14.00/feature_tools.jar

RUN curl \
        -LO \
        https://github.com/sa501428/mixer-tools/releases/download/v4.08.02/MixerTools.4.8.2.jar && \
    chmod 666 /opt/MixerTools.4.8.2.jar && \
    ln -s /opt/MixerTools.4.8.2.jar /opt/MixerTools.jar

RUN curl \
        -LO \
        https://github.com/sa501428/merge-stats/releases/download/v0.3/merge-stats.jar

RUN curl \
        -LO \
        https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip && \
    unzip -qq gatk-4.2.2.0.zip && \
    cd gatk-4.2.2.0 && \
    rm \
        -rf \
        README.md \
        gatk-completion.sh \
        gatk-package-4.2.2.0-spark.jar \
        gatkPythonPackageArchive.zip \
        gatkcondaenv.yml \
        gatkdoc \
        scripts

RUN git clone https://github.com/aidenlab/hic2gatk.git && \
    cd hic2gatk && \
    git checkout dee2a9f6f2b0e95f0cb3e8d47eef9798cb8101aa

RUN git clone --branch phasing https://github.com/aidenlab/3d-dna.git && \
    cd 3d-dna && \
    git checkout 63029aa3bc5ba9bbdad9dd9771ace583cc95e273

RUN git clone https://github.com/sa501428/psf-to-bedpe.git && \
    cd psf-to-bedpe && \
    git checkout 0.1

# For sorting, LC_ALL is C
ENV LC_ALL C
ENV PATH=/opt:/opt/scripts:/opt/scripts/common:/opt/juicer/misc:/opt/gatk-4.2.2.0:/opt/hic2gatk:$PATH

RUN mkdir -p hic-pipeline/hic_pipeline
COPY hic_pipeline hic-pipeline/hic_pipeline/
ENV PATH="/opt/hic-pipeline/hic_pipeline:${PATH}"
