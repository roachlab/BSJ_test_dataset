# Reproducible image for the BSJ benchmark:
#   STAR, BWA, samtools, CIRCexplorer2, CIRI2, Snakemake 8 (Python 3.11).
#
# Build:  docker build -t bsj-bench .
# Run:    docker run --rm -v "$PWD/results:/work/benchmark/results" bsj-bench
#         (override cores)  docker run --rm bsj-bench bash -lc "cd benchmark && snakemake --cores 8"
FROM mambaorg/micromamba:1.5.10

# 1) conda / bioconda toolchain (incl. Snakemake 8)
COPY --chown=$MAMBA_USER:$MAMBA_USER benchmark/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

# 2) CIRI2 (distributed as a Perl script, not on bioconda) -- fetched at build time.
#    Override with: --build-arg CIRI2_URL=<url>, or mount CIRI2.pl and set $CIRI2 at run time.
USER root
ARG CIRI2_URL="https://sourceforge.net/projects/ciri/files/CIRI2/CIRI_v2.0.6.zip/download"
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl unzip \
    && rm -rf /var/lib/apt/lists/* \
    && curl -L "$CIRI2_URL" -o /tmp/ciri.zip \
    && unzip -o /tmp/ciri.zip -d /opt/ciri \
    && ln -sf "$(find /opt/ciri -name 'CIRI2.pl' | head -1)" /usr/local/bin/CIRI2.pl \
    && rm -f /tmp/ciri.zip
ENV CIRI2=/usr/local/bin/CIRI2.pl

# 3) project sources
WORKDIR /work
COPY --chown=$MAMBA_USER:$MAMBA_USER . /work
USER $MAMBA_USER

# default: run the full benchmark (the micromamba entrypoint activates the env)
CMD ["bash", "-lc", "cd benchmark && snakemake --cores 4"]
