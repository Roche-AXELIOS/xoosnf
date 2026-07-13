FROM debian:13.4-slim

RUN apt-get -qq update && apt-get -qq upgrade -y \
  && apt-get -qq install -y --no-install-recommends procps python3 \
  && rm -rf /var/lib/apt/lists/*

ENV PATH=/apps:${PATH}
COPY generate-metrics-report.py /apps/generate-metrics-report
