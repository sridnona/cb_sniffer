# Published as mgibio/cb_sniffer. Build for amd64 (also from an Apple Silicon Mac):
#   docker build --platform linux/amd64 -t mgibio/cb_sniffer:<version> .
FROM python:3.11-slim

LABEL org.opencontainers.image.title="cb_sniffer" \
      org.opencontainers.image.description="Counts ref and alt cell barcodes and UMIs at variant sites in 10x single-cell BAM/CRAM files" \
      org.opencontainers.image.source="https://github.com/chrisamiller/cb_sniffer"

# procps provides `ps`, which some HPC/LSF docker wrappers expect
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /opt/cb_sniffer/requirements.txt
RUN pip install --no-cache-dir -r /opt/cb_sniffer/requirements.txt

COPY cb_sniffer.py /opt/cb_sniffer/cb_sniffer.py
RUN printf '#!/bin/sh\nexec python3 /opt/cb_sniffer/cb_sniffer.py "$@"\n' > /usr/local/bin/cb_sniffer \
    && chmod +x /usr/local/bin/cb_sniffer

# No ENTRYPOINT so the image also works with `bsub -a "docker(...)" /bin/bash`
CMD ["cb_sniffer", "-h"]
