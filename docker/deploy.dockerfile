FROM BASEDOCKERNAME
WORKDIR /opt
ENV PATH /opt/bin:$PATH
ADD data/ /opt/data/
ADD bin/ /opt/bin/
RUN diamond makedb --in /opt/data/hla.faa --db /opt/data/hla.dmnd
