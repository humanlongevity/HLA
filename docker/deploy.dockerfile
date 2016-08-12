FROM BASEDOCKERNAME
WORKDIR /opt
ENV PATH /opt/bin:$PATH
ADD data/ /opt/data/
ADD bin/ /opt/bin/
