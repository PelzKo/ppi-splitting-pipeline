# One image for every process in the pipeline, matching the single-conda-env
# style of `-profile conda` (there are no per-process containers).
#
# Build and push:
#   docker build -t docker.io/konstantinpelz/ppi-splitting:<tag> .
#   docker push  docker.io/konstantinpelz/ppi-splitting:<tag>
#
# <tag> must equal the git tag of the commit it is built from, and
# nextflow.config's `docker` profile must name the same tag. A pipeline that
# embeds this one as a git submodule pins the submodule to that tag, so the
# submodule tag and the image tag move together -- the image carries bin/ (see
# below), which means a bin/ change with a stale image is a silently wrong run.
FROM mambaorg/micromamba:2.0.5

USER root

# procps: Nextflow's trace/`ps` polling. git: cvxpy/pip metadata on some deps.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# The same environment.yml the conda profile uses, so the two profiles cannot
# drift apart in what they install. torch's cu128 pin lives in environment.yml
# itself (see the comment there) -- this is the only torch install in the
# image now; a prior version force-reinstalled a second copy here, which left
# both the unpinned and pinned wheels in the image (~18GB combined, since
# Docker layers are additive and never delete what an earlier layer wrote).
# Keep the pin in sync with the driver on the target cluster.
ENV PIP_NO_CACHE_DIR=1
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH

# Nextflow puts only the *root* project's bin/ on PATH. Standalone that is this
# repo, but when this pipeline is included as a subworkflow the root project is
# the embedding one, so sample_negatives.py et al. would not resolve. Baking them
# in makes every process work identically in both cases -- and is why the image
# tag is tied to the commit.
COPY bin/ /opt/ppi-splitting/bin/
RUN chmod -R a+rx /opt/ppi-splitting/bin
ENV PATH=/opt/ppi-splitting/bin:$PATH

USER $MAMBA_USER
