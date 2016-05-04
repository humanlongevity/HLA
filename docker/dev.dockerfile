FROM ubuntu:14.04

RUN sudo sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list' \
	&& gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 \
	&& gpg -a --export E084DAB9 | sudo apt-key add -

# [1] Get noninteractive frontend for Debian to avoid some problems:
#    debconf: unable to initialize frontend: Dialog
# [2] Always combine RUN apt-get update with apt-get install in the same RUN statement for cache busting
# https://docs.docker.com/engine/articles/dockerfile_best-practices/
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections \
 	&& sudo apt-get update -y \
 	&& sudo apt-get install -y \
		gcc \
		g++ \
		zlib1g-dev \
		libcurl4-openssl-dev \
		r-base \ 
		wget \
		python-pip \
 	&& apt-get clean \
 	&& rm -rf /var/lib/apt/lists/*		

RUN pip install --upgrade pip \
	&&  pip install awscli requests pytest boto3 logger

RUN wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 \
	&& tar xjf samtools-1.3.tar.bz2 \
	&& cd samtools-1.3 \
	&& mv htslib-1.3 temp \
	&& wget https://github.com/humanlongevity/htslib/archive/jpiper/1.3-iam-support.zip \
	&& unzip 1.3-iam-support.zip \
	&& mv htslib-jpiper-1.3-iam-support htslib-1.3 \
	&& cd htslib-1.3 && autoconf && ./configure --enable-libcurl && make -j 8 all && cd .. \
	&& perl -pi -e 's/^(LIBS\s+=)/\1 -lcurl -lcrypto/' Makefile \
	&& make -j 8 \
	&& make install \
	&& cd .. \
	&& rm -rf samtools-1.3*

RUN wget https://github.com/bbuchfink/diamond/releases/download/v0.7.12/diamond-linux64.tar.gz \
	&& tar xzf diamond-linux64.tar.gz \
	&& mv diamond /usr/bin/ \
	&& rm diamond-linux64.tar.gz

RUN echo 'install.packages("devtools", repos="http://cran.r-project.org", clean=TRUE);q()' | sudo R --no-save \
	&& echo 'devtools::install_version("data.table", "1.9.6", repos="http://cran.r-project.org", clean=TRUE);q()' | sudo R --no-save \
	&& echo 'devtools::install_version("lpSolve", "5.6.13", repos="http://cran.r-project.org", clean=TRUE);q()' | sudo R --no-save

ENTRYPOINT ["python", "/opt/bin/run.py"]
