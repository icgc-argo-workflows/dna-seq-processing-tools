FROM quay.io/icgc-argo/dna-seq-processing-tools:base-docker.0.1.0

# cramtools installation
RUN curl -sSL -o /tools/cramtools.jar --retry 10 https://github.com/enasequence/cramtools/blob/v3.0/cramtools-3.0.jar?raw=true

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

CMD ["/bin/bash"]
