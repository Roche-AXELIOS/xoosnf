FROM python:3.13-slim-trixie

RUN apt-get -qq update && apt-get -qq upgrade -y \
  && apt-get -qq install -y --no-install-recommends procps \
  && rm -rf /var/lib/apt/lists/*

ENV PATH=/apps:${PATH}
COPY subsample-run.py /apps/subsample-run
