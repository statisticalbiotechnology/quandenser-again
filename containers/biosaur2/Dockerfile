# biosaur2 MS1 feature detection, as an alternative to Dinosaur.
#
# Pinned to the upstream commit that carries the exact fast path for timsTOF
# MS1 processing (markmipt/biosaur2#22): byte-identical features, about 50x
# faster on PASEF runs. Its Cython kernels compile at install time, hence the
# compiler in the image.
#
# The fast path reproduces numpy 1.x floating-point results and falls back to
# the original, much slower code on numpy >= 2 (biosaur2/cfast.pyx, NUMPY_OK).
# biosaur2's own requirements pin numpy 1.26.4, so the image must not upgrade
# it; the build checks that the fast path is actually live.
#
#   docker build -t biosaur2:2e7c9a1 containers/biosaur2
FROM python:3.12-slim

ARG BIOSAUR2_COMMIT=2e7c9a1

RUN apt-get update && apt-get install -y --no-install-recommends \
        g++ git procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir "numpy==1.26.4" Cython setuptools wheel \
    && pip install --no-cache-dir --no-build-isolation \
        "git+https://github.com/markmipt/biosaur2.git@${BIOSAUR2_COMMIT}" \
    && python -c "import biosaur2.fastpath as f, sys; sys.exit(0 if f.FAST_OK else 'biosaur2 fast path is not active')" \
    && biosaur2 -h > /dev/null

ENTRYPOINT []
CMD ["biosaur2", "-h"]
