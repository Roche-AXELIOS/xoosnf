ARG SAMTOOLS_VERSION=1.23.1
ARG KMC_VERSION=3.2.4
ARG SEQKIT_VERSION=2.13.0

FROM cgr.dev/chainguard/wolfi-base AS kmc-builder

ARG KMC_VERSION

RUN apk add --no-cache ca-certificates curl \
  && curl -sL https://github.com/refresh-bio/KMC/releases/download/v${KMC_VERSION}/KMC${KMC_VERSION}.linux.x64.tar.gz -o /tmp/kmc.tar.gz \
  && echo "158f2084f8d928b3f33b8aaf7d1220fee4183bf46837787e5e6b16bbdf54d31d  /tmp/kmc.tar.gz" | sha256sum -c - \
  && tar zxf /tmp/kmc.tar.gz \
  && rm /tmp/kmc.tar.gz

FROM cgr.dev/chainguard/wolfi-base AS samtools-builder

ARG SAMTOOLS_VERSION

RUN apk add --no-cache build-base ca-certificates curl bzip2-dev xz-dev zlib-dev zlib-static \
  && curl -sL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 -o /tmp/samtools.tar.bz2 \
  && echo "32266198a4bc6a6df395d8526688c9697d9c8e472f888c749fdde2e08ea88dd2  /tmp/samtools.tar.bz2" | sha256sum -c - \
  && tar jxf /tmp/samtools.tar.bz2 \
  && rm /tmp/samtools.tar.bz2 \
  && cd samtools-${SAMTOOLS_VERSION} \
  && ./configure --without-curses --disable-libcurl --disable-plugins LDFLAGS="-static" \
  && make

FROM cgr.dev/chainguard/wolfi-base AS seqkit-builder

ARG SEQKIT_VERSION

RUN apk add --no-cache ca-certificates curl \
  && curl -sL https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -o /tmp/seqkit.tar.gz \
  && echo "7d686de448464fada1b1988e2e07d693bec68768312da62846bc0e2b502bfc46  /tmp/seqkit.tar.gz" | sha256sum -c - \
  && tar zxf /tmp/seqkit.tar.gz \
  && rm /tmp/seqkit.tar.gz

FROM cgr.dev/chainguard/wolfi-base

ARG SAMTOOLS_VERSION

RUN apk add --no-cache bash pigz procps python-3.13 zstd

ENV PATH="/apps:${PATH}"

WORKDIR /apps

COPY --from=kmc-builder /bin/kmc /apps/kmc
COPY --from=samtools-builder /samtools-${SAMTOOLS_VERSION}/samtools /apps/samtools
COPY --from=seqkit-builder /seqkit /apps/seqkit
COPY fastq_downsample.py /apps/fastq_downsample
