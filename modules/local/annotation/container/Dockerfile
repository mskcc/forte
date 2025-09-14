####TARGET Fusion Annotator - Kofi Amoah - 07/24/2025

################## Base Image ###########
ARG PYTHON_VERSION="3.12"
FROM python:${PYTHON_VERSION}

LABEL org.opencontainers.image.vendor="MSKCC" \
    org.opencontainers.image.authors="Kofi Amoah (amoahk1@mskcc.org)" \
    org.opencontainers.image.vcs-url="https://github.mskcc.org/ClinBX/TARGET-fusion-annotator.git" \
    org.opencontainers.image.description="This container uses python 3.12 as the base image to build \
    iAnnotateSV lite version"

################## INSTALL ##########################
WORKDIR /app
ADD . /app

RUN cd iAnnotateSV_lite \
    && apt-get update && apt-get install --no-install-recommends -y gcc g++ zlib1g-dev perl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && make deps-install \
    && poetry build \
    && pip install -r requirements.txt

ENTRYPOINT []