FROM quay.io/pancancer/dna-seq-processing:0.0.1

USER root

COPY ./tools/*/*.py /tools/

COPY cwltool-launcher/cwltool_launcher.py /tools/

CMD /bin/bash
