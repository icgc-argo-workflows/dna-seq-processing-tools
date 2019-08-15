FROM quay.io/pancancer/dna-seq-processing:0.0.1

USER root

RUN pip3 install requests

COPY ./tools/*/*.py /tools/

CMD /bin/bash
