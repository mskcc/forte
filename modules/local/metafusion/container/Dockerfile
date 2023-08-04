FROM mapostolides/metafusion

LABEL author="Alexandria Dymun (pintoa1@mskcc.org)" \
    maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.image="0.0.5"

ENV METAFUSION_TAG="v1.0.1"
ENV PATH="${PATH}:/MetaFusion/scripts"

RUN git clone -b `echo $METAFUSION_TAG` https://github.com/pintoa1-mskcc/MetaFusion.git && \
    cd /MetaFusion/ && \
    git pull
RUN chmod -R +x /MetaFusion/scripts
