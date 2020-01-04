FROM quay.io/icgc-argo/dna-seq-processing-tools:base-docker.0.1.1

# cramtools installation
# the current master has commit hash 0b5c9ece97dba616226405bec4a053e5b04c53d1
# using hash to protect us from potential future problems
RUN curl -sSL -o /tools/cramtools.jar --retry 10 https://github.com/enasequence/cramtools/blob/0b5c9ece97dba616226405bec4a053e5b04c53d1/cramtools-3.0.jar?raw=true

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

CMD ["/bin/bash"]


