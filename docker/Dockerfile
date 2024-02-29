FROM debian:bullseye-slim

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    cmake \
    libopenblas-base \
    python3-pip \
    r-base-core \
    time && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ADD ./extdata/install_packages.R /usr/local/bin/

RUN pip3 install cget && \
    Rscript /usr/local/bin/install_packages.R

WORKDIR /app

COPY . .

# Force step_2 to use 1 single thread. More threads are ineffective
ENV OMP_NUM_THREADS=1

RUN R CMD INSTALL .

RUN mv \
    extdata/step1_fitNULLGLMM_qtl.R \
    extdata/step2_tests_qtl.R \
    extdata/step3_gene_pvalue_qtl.R \
    extdata/makeGroupFile.R \
    /usr/local/bin/ && \
    chmod a+x /usr/local/bin/*.R \
    step1_fitNULLGLMM_qtl.R --help && \
    step2_tests_qtl.R --help && \
    step3_gene_pvalue_qtl.R --help && \
    makeGroupFile.R  --help
