FROM BASEDOCKERNAME
WORKDIR /opt
ENV PATH /opt/bin:$PATH
ENTRYPOINT ["python", "/opt/bin/run.py"]
ADD bin/ /opt/bin/
ADD data/ /opt/data/
