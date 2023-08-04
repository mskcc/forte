FROM mapostolides/metafusion
ENV METAFUSION_TAG="v1.0.1"
ENV PATH="${PATH}:/MetaFusion/scripts"

RUN git clone `echo $METAFUSION_TAG` https://github.com/pintoa1-mskcc/MetaFusion.git && \
    cd /MetaFusion/ && \
    git pull
RUN chmod +x /MetaFusion/scripts
