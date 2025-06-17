FROM rocker/tidyverse:4.4.3

# Install some dependencies
RUN apt-get update 
RUN apt-get install -y --no-install-recommends \
    libfftw3-dev \
    libtiff-dev \
    r-cran-irkernel

# Give the id=1000 user sudo access without password
ARG USERNAME=rstudio
RUN apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# [Optional] Set the default user. Omit if you want to keep the default as root.
USER $USERNAME

# # Install some R packages
RUN Rscript -e "install.packages(c('languageserver','devtools','IRkernel','rlist','pheatmap','feather','BiocManager'), repos='https://cloud.r-project.org', dependencies = TRUE)" 
RUN Rscript -e "BiocManager::install(c('Cardinal','CardinalWorkflows'), dependencies = TRUE)"