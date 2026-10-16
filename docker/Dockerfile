FROM mambaorg/micromamba:2.0.5
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/nfvs.yml /tmp/nfvs.yml
RUN micromamba install -y -n base -f /tmp/nfvs.yml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
