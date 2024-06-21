FROM rocker/r-ver:4.4.1

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install -y \
    build-essential \
    cmake \
    curl \
    libcurl4-openssl-dev \
    libopenblas-base \
    pigz \
    python3-pip \
    time \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install cget pandas && \
    R --slave --vanilla -e 'install.packages("pak", repos = "https://cloud.r-project.org")'

WORKDIR /app

# Force step_2 to use 1 single thread. More threads are ineffective
ENV OMP_NUM_THREADS=1

COPY DESCRIPTION .
RUN R --slave --vanilla -e 'pak::local_install_deps(dependencies = TRUE); pak::cache_clean()'

COPY . .
RUN R --slave --vanilla -e 'pak::local_install()'

RUN mv extdata/step1_fitNULLGLMM_qtl.R extdata/step2_tests_qtl.R extdata/step3_gene_pvalue_qtl.R extdata/makeGroupFile.R /usr/local/bin/ && \
    chmod a+x \
    /usr/local/bin/step1_fitNULLGLMM_qtl.R \
    /usr/local/bin/step2_tests_qtl.R \
    /usr/local/bin/step3_gene_pvalue_qtl.R \
    /usr/local/bin/makeGroupFile.R && \
    makeGroupFile.R  --help && \
    step1_fitNULLGLMM_qtl.R --help && \
    step2_tests_qtl.R --help && \
    step3_gene_pvalue_qtl.R --help

CMD ["/bin/bash"]
