FROM registry.gitlab.com/milothepsychic/r-alpine/tidyverse
LABEL maintainer="miles-smith@omrf.org"

RUN apk --no-cache --update-cache add \
	git \
	freetype-dev \
	glu-dev \
	libx11-dev \
	libxml2-dev \
	xz-dev \
	xz-libs

RUN Rscript -e "install.packages(c('pheatmap', 'kableExtra', 'IRdisplay', 'RColorBrewer', 'paletteer', 'cowplot', 'irlba', 'ggforce', 'uwot', 'formattable', 'drake', 'rsvd', 'BiocManager'))"
RUN Rscript -e "BiocManager::install(c('BiocParallel', 'clusterProfiler', 'DESeq2', 'S4Vectors', 'SummarizedExperiment', 'tximport', 'genefilter', 'gtools', 'inspectdf', 'rtracklayer', 'scales', 'utils', 'viridis'))"
RUN Rscript -e "BiocManager::install(c('atlanhq/flyio@cea5f645ce66bc5fb60ac7f7f2840fc8593da02e', 'milescsmith/moduleScoreR@e6db70fd27e3da25c1d88f839183c90d383a55e0'))"
RUN Rscript -e "BiocManager::install('R.utils')"

CMD ["Rscript"]