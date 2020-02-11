FROM milescsmith/mro:mro-rnaseq
LABEL maintainer="miles-smith@omrf.org"
RUN \
  export DEBIAN_FRONTEND=noninteractive && \
  apt-get -y update && \
  apt-get install -y \
    curl \
    git-core \
    libglu1-mesa-dev \
    pandoc \
    pandoc-citeproc \
    libcgal-dev \
    libfreetype6-dev \
    libhdf5-dev \
    libopenblas-dev \
    liblapack-dev \
    libssh2-1-dev \
    rclone \
    ssh \
    vim && \
  apt-get clean -y && \
  apt-get autoremove -y && \
  rm -rf /tmp/downloaded_packages/* && \
  rm -rf /var/lib/apt/lists/*

COPY packages.R /opt/packages.R
RUN Rscript /opt/packages.R

WORKDIR /root

CMD ["Rscript"]