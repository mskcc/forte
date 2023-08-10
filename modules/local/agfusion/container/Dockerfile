FROM ubuntu:bionic-20230530

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.image="0.0.5"

# INSTALL DEPENDENCIES

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y
RUN apt-get install -y build-essential python3 python3-pip python3-matplotlib python3-pandas python3-future python3-biopython curl less vim libnss-sss git zip
RUN pip3 install --upgrade pip
RUN pip3 install pyensembl

# Additional libraries needed for AGFusion build command
RUN apt-get install -y default-libmysqlclient-dev
RUN pip3 install mysqlclient

# INSTALL AGFUSION & DATABASE FILES
WORKDIR /usr/local/bin
RUN git clone https://github.com/anoronh4/AGFusion.git --branch feature/metafusion_parser
WORKDIR /usr/local/bin/AGFusion
RUN pip3 install .

# downgrade pyensembl for compatibility
RUN pip3 install gtfparse==1.2.1 --upgrade
