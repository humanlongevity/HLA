NS = docker-dev.hli.io/xchao
TAG = 0.6.1
NAME = hla

build:
	docker build -t $(NS)/$(NAME)-dev -f dev.dockerfile .
	docker tag $(NS)/$(NAME)-dev $(NS)/$(NAME)-dev:$(TAG)

deploy:
	sed "s|BASEDOCKERNAME|$(NS)/$(NAME)-dev:$(TAG)|" < deploy.dockerfile > deploy.filled.dockerfile
	rsync -avP ../bin ./
	rsync -avP ../data/hla.* ./data/
	rsync -avP ../data/chr6 ./data/
	docker build -t $(NS)/$(NAME) -f deploy.filled.dockerfile .
	docker tag $(NS)/$(NAME) $(NS)/$(NAME):$(TAG)
	rm -r bin/ data/

push: build deploy
	docker push $(NS)/$(NAME)-dev
	docker push $(NS)/$(NAME)-dev:$(TAG)
	docker push $(NS)/$(NAME)
	docker push $(NS)/$(NAME):$(TAG)

shell:
	docker run --rm -ti \
		-v ${HOME}/.aws:/root/.aws:ro \
		-v ${TMPDIR}:/scratch \
		-v `pwd`:/work \
		-w /work \
		-e TMPDIR=/scratch \
		--net=host \
		--entrypoint='/bin/bash' \
		$(NS)/$(NAME)-dev:$(TAG)

shell2:
	docker run --rm -ti \
		-v ${HOME}/.aws:/root/.aws:ro \
		-v ${TMPDIR}:/scratch \
		-v `pwd`:/work \
		-w /work \
		-e TMPDIR=/scratch \
		--net=host \
		--entrypoint='/bin/bash' \
		$(NS)/$(NAME):$(TAG)

test-integration:
	python integration_tests.py --bind-source --image=$(NS)/$(NAME)-dev:$(TAG)

test-deploy:
	python integration_tests.py --image=$(NS)/$(NAME):$(TAG)

clean:
	rm -rf _candidates deploy.filled.dockerfile
