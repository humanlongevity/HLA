NS = docker-dev.hli.io/xchao
VERSION = 0.3

NAME = hla
BAM_FILE = s3://hli-processed/processed/220033/isis_alignment_and_variant_calling/176444255/Data/Intensities/BaseCalls/Alignment/176444255_S1.bam
SAMPLE_ID = 176444255

.PHONY: all
all: build push

build:
	rsync -avP ../bin ../data .
	docker build -t $(NS)/$(NAME):$(VERSION) -f Dockerfile . 

push:
	docker push $(NS)/$(NAME):$(VERSION)

test:
	./hla.sh $(NS)/$(NAME):$(VERSION) $(BAM_FILE) $(SAMPLE_ID)
