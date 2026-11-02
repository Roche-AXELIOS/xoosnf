FROM cgr.dev/chainguard/wolfi-base

RUN apk add --no-cache bash python-3.13 procps

ENV PATH=/apps:${PATH}
COPY parsers.py /apps/parsers.py
COPY xoos_to_multiqc.py /apps/xoos_to_multiqc.py
COPY combine_reports.py /apps/combine_reports.py
