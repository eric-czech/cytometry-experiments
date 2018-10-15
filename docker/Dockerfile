FROM gfinak/opencyto:latest


######################
# Shiny Installation #
######################
#
# From: https://github.com/rocker-org/shiny/blob/master/Dockerfile

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget

# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('shiny', 'rmarkdown'), repos='$MRAN')" && \
    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/

# Install bioconductor tools
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("FlowRepositoryR", update=F)'
# RUN R -e 'BiocManager::install("flowDensity", update=F)'

# Install packages for shinyCyto
RUN R -e 'devtools::install_github("hammerlab/cytoEx", ref="trunk")'
RUN R -e 'devtools::install_github("RGLab/flowWorkspaceData", ref="trunk")'
RUN R -e "devtools::install_github('AnalytixWare/ShinySky')"
RUN R -e "install.packages(c('shinyjs', 'shinydashboard', 'DT', 'shinyFiles', 'networkD3', 'shinyBS', 'hash', 'snakecase'))"

EXPOSE 3838

COPY shiny-server.sh /usr/bin/shiny-server.sh