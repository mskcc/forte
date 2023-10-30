FROM ubuntu:14.04

LABEL maintainer="Anne Marie (noronhaa@mskcc.org)" \
    version.image="1.2.4" \
    source.getBaseCountsMultiSample="https://github.com/zengzheng123/GetBaseCountsMultiSample/releases/tag/v1.2.2" \
    version.getBaseCountsMultiSample="1.2.4"
ENV GBCMS_VERSION 1.2.4

RUN apt-get update && \
    apt-get install --yes \
    procps \
    gcc \
    make \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    g++ \
    git \
    wget \
    zip \
    python-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN cd /tmp && \
    wget https://github.com/zengzheng123/GetBaseCountsMultiSample/archive/v${GBCMS_VERSION}.zip && \
    unzip v${GBCMS_VERSION}.zip && \
    # Install bamtools
    cd /tmp/GetBaseCountsMultiSample-${GBCMS_VERSION}/bamtools-master && \
    rm -r build/ && \
    mkdir build && \
    cd build/ && \
    cmake -DCMAKE_CXX_FLAGS=-std=c++03 .. && \
    make && \
    make install && \
    cp ../lib/libbamtools.so.2.3.0 /usr/lib/ && \
    # Install GetBaseCountsMultiSample itself
    cd /tmp/GetBaseCountsMultiSample-${GBCMS_VERSION} && \
    make && \
    cp GetBaseCountsMultiSample /usr/bin/
