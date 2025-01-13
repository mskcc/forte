FROM ubuntu:jammy-20240911.1

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.image="0.0.7"

# INSTALL DEPENDENCIES

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y
RUN apt-get install -y \
    build-essential \
    python3 \
    python3-pip \
    python3-matplotlib \
    python3-pandas \
    python3-future \
    python3-biopython \
    python3-dev \
    default-libmysqlclient-dev \
    pkg-config \
    curl \
    less \
    vim \
    libnss-sss \
    git \
    zip
RUN pip3 install --upgrade pip
RUN pip3 install pyensembl

# Additional libraries needed for AGFusion build command
RUN apt-get install -y default-libmysqlclient-dev
RUN pip3 install mysqlclient

# INSTALL AGFUSION & DATABASE FILES
WORKDIR /usr/local/bin
RUN git clone https://github.com/mskcc/AGFusion.git --branch v1.4.3@mskcc.1 --single-branch
WORKDIR /usr/local/bin/AGFusion
RUN pip3 install -r requirements.txt
RUN pip3 install .

