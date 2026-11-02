FROM cgr.dev/chainguard/wolfi-base

RUN apk add --no-cache bash python-3.13 procps

ENV PATH=/apps:${PATH}
COPY generate-metrics-report.py /apps/generate-metrics-report
