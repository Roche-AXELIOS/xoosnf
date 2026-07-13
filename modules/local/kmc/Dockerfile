ARG SAMTOOLS_VERSION=1.23.1
ARG KMC_VERSION=3.2.4

FROM debian:13.4-slim AS kmc-builder

ARG KMC_VERSION

RUN apt-get -qq update && apt-get -qq install -y --no-install-recommends ca-certificates curl \
  && rm -rf /var/lib/apt/lists/* \
  && curl -sL https://github.com/refresh-bio/KMC/releases/download/v${KMC_VERSION}/KMC${KMC_VERSION}.linux.x64.tar.gz -o /tmp/kmc.tar.gz \
  && echo "158f2084f8d928b3f33b8aaf7d1220fee4183bf46837787e5e6b16bbdf54d31d  /tmp/kmc.tar.gz" | sha256sum -c - \
  && tar zxf /tmp/kmc.tar.gz \
  && rm /tmp/kmc.tar.gz

FROM debian:13.4-slim AS samtools-builder

ARG SAMTOOLS_VERSION

RUN apt-get -qq update && apt-get -qq install -y --no-install-recommends \
    build-essential ca-certificates curl libbz2-dev liblzma-dev libncurses-dev zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && curl -sL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 -o /tmp/samtools.tar.bz2 \
  && echo "32266198a4bc6a6df395d8526688c9697d9c8e472f888c749fdde2e08ea88dd2  /tmp/samtools.tar.bz2" | sha256sum -c - \
  && tar jxf /tmp/samtools.tar.bz2 \
  && rm /tmp/samtools.tar.bz2 \
  && cd samtools-${SAMTOOLS_VERSION} \
  && ./configure --without-curses --disable-libcurl --disable-plugins LDFLAGS="-static" \
  && make

FROM staphb/seqkit:2.13.0 AS seqkit-builder

FROM debian:13.4-slim

ARG SAMTOOLS_VERSION

RUN apt-get -qq update && apt-get -qq upgrade -y \
  && apt-get -qq install -y --no-install-recommends \
    pigz \
    procps \
    python3 \
    zstd \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="/apps:${PATH}"

WORKDIR /apps

COPY --from=kmc-builder /bin/kmc /apps/kmc
COPY --from=samtools-builder /samtools-${SAMTOOLS_VERSION}/samtools /apps/samtools
COPY --from=seqkit-builder /usr/local/bin/seqkit /apps/seqkit
COPY fastq_downsample.py /apps/fastq_downsample
