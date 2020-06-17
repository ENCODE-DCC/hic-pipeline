FROM ubuntu:16.04@sha256:e9938f45e51d9ff46e2b05a62e0546d0f07489b7f22fbc5288defe760599e38a as main

LABEL maintainer "Paul Sud"
LABEL maintainer.email "encode-help@lists.stanford.edu"

# Package versions
ENV BWA_VERSION=0.7.17
ENV JUICER_TOOLS_VERSION=1.11.09_jcuda.0.8
ENV PAIRIX_VERSION=0.3.6
ENV SAMTOOLS_VERSION=1.6

RUN apt-get update && apt-get install -y \
    bzip2 \
    curl \
    gawk \
    gcc \
    git \
    libbz2-dev \
    libz-dev \
    locales \
    make \
    openjdk-8-jdk \
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

# Install pairix, only need bam2pairs + deps
RUN curl -OL "https://github.com/4dn-dcic/pairix/archive/${PAIRIX_VERSION}.zip" && \
    unzip "${PAIRIX_VERSION}.zip" && \
    cd "pairix-${PAIRIX_VERSION}/" && \
    make && \
    cp bin/pairix bin/bgzip /usr/local/bin && \
    cp util/bam2pairs/bam2pairs /opt && \
    cd .. && \
    rm -rf "${PAIRIX_VERSION}.zip" "pairix-${PAIRIX_VERSION}"

# Install Juicer
RUN git clone --branch encode https://github.com/theaidenlab/juicer.git && \
    cd juicer && \
    git checkout da9b33227756f8b00d831fc52f951baa984a68aa && \
    chmod +x CPU/* CPU/common/* misc/* && \
    find -mindepth 1 -maxdepth 1  -type d -not -name "CPU" -not -name ".git" -not -name "misc" | xargs rm -rf

# Install Juicer tools
RUN curl -L "http://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_${JUICER_TOOLS_VERSION}.jar" -o /opt/juicer/CPU/common/juicer_tools.jar && \
    chmod 666 /opt/juicer/CPU/common/juicer_tools.jar && \
    ln -s juicer/CPU scripts

# For sorting, LC_ALL is C
ENV LC_ALL C
ENV PATH=/opt:/opt/scripts:/opt/scripts/common:/opt/juicer/misc:$PATH

# Use multistage build for fragment
FROM rust:1.42.0-slim-stretch@sha256:1201e100253c82c4accd19c649d187ab741b89954397dcdddb792008f898f93a as builder
WORKDIR /hic-pipeline
COPY . .
RUN cargo build --release

FROM main

RUN mkdir -p hic-pipeline/hic_pipeline
COPY hic_pipeline hic-pipeline/hic_pipeline/
RUN rm -rf hic-pipeline/hic_pipeline/fragment
# Add compiled binary from other stage of build
COPY --from=builder /hic-pipeline/target/release/fragment hic-pipeline/hic_pipeline
ENV PATH="/opt/hic-pipeline/hic_pipeline:${PATH}"
