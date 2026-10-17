FROM --platform=linux/amd64 python:3.11-slim

# procps provides `ps`, which some HPC/LSF docker wrappers expect
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

# pysam 0.14.1 (requirements.txt) has no wheels for modern Python; the APIs
# used by cb_sniffer are unchanged in current releases
RUN pip install --no-cache-dir pysam==0.22.1

COPY cb_sniffer.py /opt/cb_sniffer/cb_sniffer.py
RUN printf '#!/bin/sh\nexec python3 /opt/cb_sniffer/cb_sniffer.py "$@"\n' > /usr/local/bin/cb_sniffer \
    && chmod +x /usr/local/bin/cb_sniffer

# No ENTRYPOINT so the image also works with `bsub -a "docker(...)" /bin/bash`
CMD ["cb_sniffer", "-h"]
