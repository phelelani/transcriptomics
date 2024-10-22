Bootstrap:shub
From:singularityhub/ubuntu

%labels
Maintainer Phelelani.Mpangase@wits.ac.za

%post
## Updates and essentials
apt-get update
apt-get install -y build-essential
apt-get install -y wget git curl gfortran zlib1g-dev libbz2-dev liblzma-dev libpcre3-dev libcurl4-gnutls-dev
apt-get install -y libxml2-dev libssl-dev libopenblas-dev libmagick++-dev libudunits2-dev libcairo2-dev libxt-dev

## Install R
cd /opt \
    && wget https://cran.r-project.org/src/base/R-3/R-3.5.1.tar.gz \
    && tar -vxf R-3.5.1.tar.gz \
    && cd R-3.5.1 \
    && ./configure --with-x=no --with-readline=no \
    && make \
    && make install \
    && rm /opt/R-3.5.1.tar.gz

## Install basic CRAN R packages
R -e 'install.packages(c("tidyverse","data.table","dtplyr","devtools","roxygen2", "Cairo"), repos="http://cloud.r-project.org/", dependencies=TRUE)'

## Install CRAN R packages
R -e 'install.packages(c("ggplot2", "gridExtra", "ggrepel", "xtable", "gplots", "kableExtra", "grid", "pheatmap", "enrichR", "BiocManager", "UpSetR", "PoiClaClu"), repos="http://cloud.r-project.org/", dependencies=TRUE)'

## Install Bioconductor R packages
R -e 'BiocManager::install(c("DESeq2", "biomaRt", "genefilter", "AnnotationDbi", "org.Hs.eg.db", "pathview", "gage", "gageData", "PROPER"), version = "3.8")'
