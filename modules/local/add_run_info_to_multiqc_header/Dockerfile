FROM cgr.dev/chainguard/wolfi-base

RUN apk add --no-cache bash python-3.13 procps

ENV PATH=/apps:${PATH}
COPY add_run_info_to_multiqc_header.py /apps/add-run-info-to-multiqc-header
