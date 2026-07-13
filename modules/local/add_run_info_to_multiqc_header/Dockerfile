FROM debian:13.4-slim

RUN apt-get -qq update && apt-get -qq upgrade -y \
  && apt-get -qq install -y --no-install-recommends procps python3 \
  && rm -rf /var/lib/apt/lists/*

ENV PATH=/apps:${PATH}
COPY add_run_info_to_multiqc_header.py /apps/add-run-info-to-multiqc-header
