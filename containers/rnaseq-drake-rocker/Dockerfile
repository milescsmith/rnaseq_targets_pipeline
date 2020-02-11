FROM bioconductor/release_core2:R3.6.1_Bioc3.10

LABEL maintainer="miles-smith@omrf.org"

RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update && \
	apt-get install -y git-core \
		tzdata \
		libbz2-dev \
		libcurl4-openssl-dev \
		liblzma-dev \
		libpng-dev \
		libssl-dev \
		libxml2-dev \
		libglu1-mesa-dev \
		libx11-dev \
		make \
		pandoc \
		pandoc-citeproc \
		python-minimal && \
	cp /usr/share/zoneinfo/America/Chicago /etc/localtime && \
	apt-get clean && \
	rm -rf /tmp/downloaded_packages/* && \
	rm -rf /var/lib/apt/lists/*

USER rstudio

COPY packages.R /opt/packages.R
RUN Rscript /opt/packages.R

WORKDIR /home/rstudio

USER root

EXPOSE 8787
#CMD ["Rscript"]