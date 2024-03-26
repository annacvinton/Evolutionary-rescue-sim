IMAGE_ROOT?=ghcr.io/erikkusch
IMAGE=evores
IMAGE_VERSION=0.2.1
R_VERSION=4.3.2


build: Dockerfile
	docker buildx build --platform linux/amd64 \
		--label "org.opencontainers.image.source=https://github.com/annacvinton/Evolutionary-rescue-sim" \
		--label "org.opencontainers.image.description=EvoRes environment with R $(R_VERSION)" \
		--build-arg R_VERSION=$(R_VERSION) \
		-t $(IMAGE_ROOT)/$(IMAGE):$(IMAGE_VERSION) \
		.

push:
	docker push $(IMAGE_ROOT)/$(IMAGE):$(IMAGE_VERSION)

singularity:
	rm -f $(IMAGE).sif $(IMAGE).tar
	docker save $(IMAGE_ROOT)/$(IMAGE):$(IMAGE_VERSION) -o $(IMAGE).tar
	singularity build $(IMAGE).sif docker-archive://$(IMAGE).tar
	rm -f $(IMAGE).tar
