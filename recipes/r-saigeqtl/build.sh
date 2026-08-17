#!/bin/bash

if [[ "$(uname -m)" == "arm64" ]]; then
    # let's dynamically detect the GCC and Darwin version in the gfortran library path
    GCC_LIB_PATH=$(find "${PREFIX}/lib/gcc/" -type d -name "arm64-apple-darwin*" | head -n 1)
    GCC_VERSION=$(basename "$(find "${GCC_LIB_PATH}" -mindepth 1 -maxdepth 1 -type d | sort | tail -n 1)")
    # hard gfortran paths in R build(?) break the build on osx-arm64
    sed -i.bak "s|\\\$(FLIBS)|-L${GCC_LIB_PATH}/${GCC_VERSION} -lgfortran -lquadmath -lm -L\${PREFIX}/lib/R/lib -lR -Wl,-framework -Wl,CoreFoundation|" src/Makevars
fi

$R CMD INSTALL --build .

pushd extdata
install "step1_fitNULLGLMM_qtl.R" "step2_tests_qtl.R" "step2b_QC_qtl.R" "step3_gene_pvalue_qtl.R" "makeGroupFile.R" "${PREFIX}/bin"
popd
