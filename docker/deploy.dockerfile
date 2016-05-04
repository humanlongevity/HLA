FROM BASEDOCKERNAME
WORKDIR /opt
ENV PATH /opt/bin:$PATH
ADD bin/ /opt/bin/
ADD data/ /opt/data/
