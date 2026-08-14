FROM mambaorg/micromamba:2.0.5

COPY --chown=$MAMBA_USER:$MAMBA_USER container/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/mashpit-isolate-screen
ENV PATH="/opt/mashpit-isolate-screen/scripts:${PATH}"
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python3", "/opt/mashpit-isolate-screen/scripts/screen_isolate.py"]
