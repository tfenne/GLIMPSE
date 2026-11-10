---
layout: default
title: Docker
nav_order: 3
parent: Installation
has_children: false
permalink: /docs/installation/docker
---

# Docker image

The `Dockerfile` at the root of the repository builds an image with all the GLIMPSE2 tools (`GLIMPSE2_chunk`, `GLIMPSE2_split_reference`, `GLIMPSE2_phase`, `GLIMPSE2_ligate`, `GLIMPSE2_concordance` and `GLIMPSE2_inspect`), plus `bcftools`, `bgzip` and `tabix`, on Ubuntu 26.04 LTS. It builds for `linux/amd64` and `linux/arm64` (e.g. AWS Graviton or Apple silicon).

The amd64 image needs a CPU with AVX2 and FMA; on Linux, `grep -ow -e avx2 -e fma /proc/cpuinfo | sort -u` should print both.

## Building the image

Clone the repository with its SIMDe submodule:

```
git clone --recursive https://github.com/odelaneau/GLIMPSE.git
cd GLIMPSE
```

(In an existing clone, run `git submodule update --init` instead.)

To build an image for the architecture of the machine you are on:

```
docker build -t glimpse:2.1.0 \
  --build-arg GLIMPSE_COMMIT=$(git rev-parse --short HEAD) \
  --build-arg GLIMPSE_DATE=$(git log -1 --format=%cd --date=short) .
```

The two build arguments set the commit and date that the tools report; the image has no git history to read them from. Without them, the tools report `unknown`.

To build both architectures as one multi-architecture image and push it to a registry:

```
docker buildx build --platform linux/amd64,linux/arm64 -t <registry>/glimpse:2.1.0 --push \
  --build-arg GLIMPSE_COMMIT=$(git rev-parse --short HEAD) \
  --build-arg GLIMPSE_DATE=$(git log -1 --format=%cd --date=short) .
```

This needs a buildx builder that can build for both platforms. Docker can emulate the other architecture (QEMU, or Rosetta in Docker Desktop on Apple silicon), but compiling under emulation is much slower: on an Apple M4, the amd64 image takes about 10 minutes under emulation, against under a minute for the native arm64 image. A builder with one native node per architecture avoids that (see `docker buildx create --append`). Without `--push`, `--load` keeps the multi-architecture image locally when Docker uses the containerd image store.

## Running the tools

```
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/data -w /data glimpse:2.1.0 \
  GLIMPSE2_phase --bam-file NA12878.bam --reference ref_chr22_16000000_20000000.bin --output out.bcf
```

The tools are on the `PATH`. `--user` makes output files belong to you rather than root.

BAM/CRAM and VCF/BCF files can be read from `https://`, `s3://` and `gs://` URLs through HTSlib. Binary reference panels (`.bin`), genetic maps and other text inputs must be local files.

## Build arguments

| Argument | Default | Meaning |
|---|---|---|
| `GLIMPSE_COMMIT`, `GLIMPSE_DATE` | `unknown` | Commit and date the tools report. |
| `UBUNTU_VERSION` | `26.04` | Ubuntu release for both the build and the runtime stages. |
| `GCC_VERSION` | `15` | g++ used on amd64. `16` selects Ubuntu 26.04's GCC 16 snapshot. |
| `LLVM_VERSION` | `23` | clang used on arm64, from [apt.llvm.org](https://apt.llvm.org). clang builds a much faster `GLIMPSE2_phase` than g++ on arm64. |
| `HTSLIB_VERSION`, `HTSLIB_SHA256` | `1.25` and its SHA-256 | HTSlib, built from source. Set both together. Avoid 1.22 and 1.23, with which `GLIMPSE2_concordance` can crash. |
| `BCFTOOLS_VERSION`, `BCFTOOLS_SHA256` | `1.25` and its SHA-256 | bcftools, built from source against the same HTSlib. Set both together. |

## Binary reference panels and Boost

Binary reference panels (`.bin`) are Boost serialization archives, and Boost refuses to read an archive written by a newer archive version. The image uses Ubuntu 26.04's Boost 1.90, which writes archive version 20, like every Boost from 1.84 onwards (including the Boost 1.85 of the current Bioconda GLIMPSE2 2.0.1). Panels written by the image can be read by those builds and the image can read theirs, but builds with Boost 1.83 or older, such as the v2.0.0 Docker image, cannot read them.
