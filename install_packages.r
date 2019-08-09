# Install R packages required for GLS analysis
pkgs = c('openssl', 'httr', 'RNeXML', 'magick', 'ade4', 'phytools', 'rredis', 'ggplot2', 'reshape', 'phytools', 'BiocManager', 'phylobase', 'phylosignal', 'dplyr', 'kde1d', 'gridExtra', 'minerva', 'lmtest', 'Metrics')
#pkgs = c('spdep', 'qrng', 'adegenet', 'cctools', 'adephylo', 'kde1d')
#pkgs = c('openssl', 'httr', 'RNeXML')
install.packages(pkgs, repos='https://cloud.r-project.org')
BiocManager::install(c('rhdf5'))
