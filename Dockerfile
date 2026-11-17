FROM mambaorg/micromamba:2.0.5

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER requirements-tools.txt /tmp/requirements-tools.txt
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y
RUN PYTHONNOUSERSITE=1 micromamba run -n base python -m pip install --upgrade pip \
    && PYTHONNOUSERSITE=1 micromamba run -n base python -m pip install -r /tmp/requirements-tools.txt

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/gwasfinemapping
WORKDIR /opt/gwasfinemapping

ENV RENV_PATHS_LIBRARY=/opt/gwasfinemapping-r-lib \
    RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
RUN mkdir -p /opt/gwasfinemapping-r-lib \
    && Rscript -e 'renv::restore(prompt = FALSE)' \
    && R_LIB="$(Rscript -e 'cat(.libPaths()[1])')" \
    && R_PROFILE_USER="" R CMD INSTALL --no-multiarch --library="$R_LIB" .

USER $MAMBA_USER
ENV PATH=/opt/conda/bin:$PATH \
    JAVA_HOME=/opt/conda/lib/jvm \
    PYTHONNOUSERSITE=1 \
    GWASFM_IN_CONTAINER=1 \
    RENV_PROJECT=/opt/gwasfinemapping \
    R_PROFILE_USER=/opt/gwasfinemapping/.Rprofile \
    RENV_PATHS_LIBRARY=/opt/gwasfinemapping-r-lib \
    RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
