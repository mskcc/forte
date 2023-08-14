FROM mapostolides/metafusion

LABEL author="Alexandria Dymun (pintoa1@mskcc.org)" \
    maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.image="0.0.6"

ENV METAFUSION_TAG="v1.0.1"
ENV PATH="${PATH}:/MetaFusion/scripts"

RUN R -e "chooseCRANmirror(ind=52); install.packages(c('plyr','data.table'))"

RUN git clone -b $METAFUSION_TAG https://github.com/pintoa1-mskcc/MetaFusion.git --single-branch && \
    cd /MetaFusion/ && \
    git pull && \
    chmod -R +x /MetaFusion/scripts && \
    rm -rf /MetaFusion/test_data

# run cleanup
RUN conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

