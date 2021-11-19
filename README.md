sagvc
=====

Somatic and Germline Variant Calling Pipeline

[![Test](https://github.com/dceoy/sagvc/actions/workflows/test.yml/badge.svg)](https://github.com/dceoy/sagvc/actions/workflows/test.yml)
[![CI to Docker Hub](https://github.com/dceoy/sagvc/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/dceoy/sagvc/actions/workflows/docker-publish.yml)

Installation
------------

```sh
$ pip install -U https://github.com/dceoy/sagvc/archive/main.tar.gz
```

Dependent commands:

- `wget`
- `pigz`
- `pbzip2`
- `bgzip`
- `tabix`
- `samtools`
- `bedtools`
- `bcftools`
- `java`
- `gatk`
- `python` (Python 3)
- `python2`
- `R`
- `Rscript`
- `configManta.py`
- `delly`
- `cnvkit.py`
- `msisensor` (or `msisensor-pro`)

Docker image
------------

Pull the image from [Docker Hub](https://hub.docker.com/r/dceoy/sagvc/).

```sh
$ docker image pull dceoy/sagvc
```
