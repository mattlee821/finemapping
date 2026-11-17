FROM mambaorg/micromamba:2.0.5

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER requirements-tools.txt /tmp/requirements-tools.txt
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y
RUN micromamba run -n base python -m pip install --upgrade pip \
    && micromamba run -n base python -m pip install -r /tmp/requirements-tools.txt

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/gwasfinemapping
WORKDIR /opt/gwasfinemapping

RUN Rscript -e 'renv::restore(prompt = FALSE)' \
    && R CMD INSTALL .

USER $MAMBA_USER
ENV PATH=/opt/conda/bin:$PATH \
    JAVA_HOME=/opt/conda/lib/jvm \
    PYTHONNOUSERSITE=1
