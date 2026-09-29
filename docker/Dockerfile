# One image for every process in the pipeline, matching the single-conda-env
# style of `-profile conda` (there are no per-process containers).
#
# Build and push -- from the repository root, not from docker/: the build context
# has to contain environment.yml and bin/, which the COPY lines below reach.
#   docker build -f docker/Dockerfile -t docker.io/konstantinpelz/ppi-splitting:<tag> .
#   docker push  docker.io/konstantinpelz/ppi-splitting:<tag>
#
# <tag> must equal the git tag of the commit it is built from, and
# nextflow.config's `docker` profile must name the same tag. A pipeline that
# embeds this one as a git submodule pins the submodule to that tag, so the
# submodule tag and the image tag move together -- the image carries bin/ (see
# below), which means a bin/ change with a stale image is a silently wrong run.
FROM mambaorg/micromamba:2.0.5

USER root

# procps: Nextflow's trace/`ps` polling. ca-certificates: HTTPS for the pip
# installs below and for the UniProt/Pfam/Gurobi WLS calls tasks make.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# environment.yml is the only package source: every conda and pip dependency,
# including torch and its CUDA build (the default PyPI wheel, CUDA 13.0 for
# the pinned version -- see the comment in that file), comes from this one install. Do not
# add packages or a separate torch/CUDA layer here -- the image must install
# exactly what `-profile conda` installs, and a second torch layer once shipped
# two wheels (~18GB, Docker layers never delete what an earlier one wrote).
# To change a version, edit environment.yml and rebuild.
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
