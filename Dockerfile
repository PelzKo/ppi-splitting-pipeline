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
# drift apart in what they install.
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH

# environment.yml:19 is a bare `torch`, so pip resolves whatever wheel is current
# at build time. That is the exact failure EMBED_SEQUENCES' --require-gpu exists to
# catch: an unpinned install once fetched a cu13x wheel onto a CUDA 12.8 driver,
# torch.cuda.is_available() returned False, and the task ran 12 h on CPU before
# SLURM killed it. The rule is that the wheel's CUDA major must not exceed the
# driver's (minors are compatible), so the image pins one wheel explicitly and
# overwrites whatever the environment.yml step resolved.
#
# Keep this in sync with the driver on the target cluster. Once environment.yml
# carries the pin itself (docs/ddi-review-plan.md, "Outstanding after Phase 3"),
# delete this layer rather than maintaining two sources of truth.
ARG TORCH_VERSION=2.10.0
ARG TORCH_CUDA=cu128
RUN pip install --no-cache-dir --force-reinstall \
      --extra-index-url https://download.pytorch.org/whl/${TORCH_CUDA} \
      torch==${TORCH_VERSION}+${TORCH_CUDA}

# Nextflow puts only the *root* project's bin/ on PATH. Standalone that is this
# repo, but when this pipeline is included as a subworkflow the root project is
# the embedding one, so sample_negatives.py et al. would not resolve. Baking them
# in makes every process work identically in both cases -- and is why the image
# tag is tied to the commit.
COPY bin/ /opt/ppi-splitting/bin/
RUN chmod -R a+rx /opt/ppi-splitting/bin
ENV PATH=/opt/ppi-splitting/bin:$PATH

USER $MAMBA_USER
