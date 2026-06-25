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
# Put the base env on PATH so tools resolve in any shell (login or not).
ENV PATH=/opt/conda/bin:$PATH

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

# 2b) Rebuild STAR from source WITHOUT AVX2 so it runs on older CPUs (e.g. Ivy Bridge / E5 v2).
#     The bioconda STAR is AVX2-compiled and SIGILLs ("illegal instruction") on such hosts.
#     SSE4.2 runs everywhere; the perf cost is negligible on this tiny genome.
#     For an even older CPU use --build-arg STAR_SIMD=sse.
ARG STAR_VERSION=2.7.11b
ARG STAR_SIMD="-msse4.2"
RUN apt-get update && apt-get install -y --no-install-recommends build-essential zlib1g-dev curl ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && curl -L "https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz" -o /tmp/star.tgz \
    && tar -xzf /tmp/star.tgz -C /tmp \
    && make -C "/tmp/STAR-${STAR_VERSION}/source" STAR CXXFLAGS_SIMD="${STAR_SIMD}" \
    && cp "/tmp/STAR-${STAR_VERSION}/source/STAR" /opt/conda/bin/STAR \
    && rm -rf /tmp/star.tgz "/tmp/STAR-${STAR_VERSION}"

# 3) project sources
WORKDIR /work
COPY --chown=$MAMBA_USER:$MAMBA_USER . /work
USER $MAMBA_USER

# default: run the full benchmark (the micromamba entrypoint activates the env)
CMD ["bash", "-c", "cd benchmark && snakemake --cores 4"]
