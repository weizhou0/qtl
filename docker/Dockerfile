FROM ubuntu:20.04

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    cmake \
    libopenblas-base \
    python3-pip \
    r-cran-devtools

RUN pip3 install cget

ADD ./docker/install_packages.R /usr/local/bin/

RUN Rscript /usr/local/bin/install_packages.R

WORKDIR /app


COPY . .

#RUN Rscript extdata/install_packages.R

# Force step_2 to use 1 single thread. More threads are ineffective
ENV OMP_NUM_THREADS=1

RUN R CMD INSTALL .

RUN mv extdata/step1_fitNULLGLMM_qtl.R extdata/step2_tests_qtl.R extdata/step3_gene_pvalue_qtl.R extdata/makeGroupFile.R /usr/local/bin/

RUN chmod a+x /usr/local/bin/step1_fitNULLGLMM_qtl.R
RUN chmod a+x /usr/local/bin/step2_tests_qtl.R
RUN chmod a+x /usr/local/bin/step3_gene_pvalue_qtl.R
RUN chmod a+x /usr/local/bin/makeGroupFile.R

RUN makeGroupFile.R  --help
RUN step1_fitNULLGLMM_qtl.R --help
RUN step2_tests_qtl.R --help
RUN step3_gene_pvalue_qtl.R --help

RUN apt-get update
RUN apt-get install time
