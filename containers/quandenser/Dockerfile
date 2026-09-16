# Quandenser container image.
#
# Two stages: a builder that compiles ProteoWizard, Boost, MaRaCluster,
# Percolator and Quandenser from source and builds the Dinosaur jar with
# Maven, and a slim runtime that carries only the binary, the jar and a JRE.
#
# Build:
#   docker build -t quandenser:0.03.2 -f containers/quandenser/Dockerfile .
#
# The build context must be the root of the quandenser repository with its
# submodules checked out (git submodule update --init --recursive).
#
# Behind a TLS-inspecting proxy, drop the CA certificate into
# containers/quandenser/certs/ and pass the usual HTTP_PROXY/HTTPS_PROXY
# build args; both are no-ops otherwise.

# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------
FROM ubuntu:22.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_JOBS=4

# Optional extra CA certificates, for builds behind a TLS-inspecting proxy.
# The directory normally holds nothing but .gitkeep, in which case every step
# below is inert and apt keeps talking to the default archives over plain HTTP.
# When a bundle is supplied we also move apt to HTTPS, since such proxies
# typically only relay CONNECT.
COPY containers/quandenser/certs/ /tmp/extra-certs/
RUN set -eu; \
    if ls /tmp/extra-certs/*.crt >/dev/null 2>&1; then \
      mkdir -p /etc/ssl/certs /usr/local/share/ca-certificates/extra; \
      cat /tmp/extra-certs/*.crt >> /etc/ssl/certs/ca-certificates.crt; \
      cp /tmp/extra-certs/*.crt /usr/local/share/ca-certificates/extra/; \
      sed -i 's|http://|https://|g' /etc/apt/sources.list; \
    fi

RUN apt-get update && apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        cmake \
        file \
        g++ \
        make \
        maven \
        openjdk-8-jdk \
        patch \
        python3 \
        rsync \
        tar \
        wget \
        xz-utils \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Teach the JDK about the same extra certificates, so Maven can reach Central.
RUN set -eu; \
    for crt in /usr/local/share/ca-certificates/extra/*.crt; do \
      [ -e "$crt" ] || continue; \
      update-ca-certificates; \
      csplit -z -f /tmp/split-ca- -b '%03d.crt' "$crt" '/-----BEGIN CERTIFICATE-----/' '{*}' >/dev/null; \
      for one in /tmp/split-ca-*.crt; do \
        keytool -importcert -noprompt -trustcacerts \
                -alias "extra-$(basename "$one" .crt)" \
                -file "$one" \
                -keystore /usr/lib/jvm/java-8-openjdk-amd64/jre/lib/security/cacerts \
                -storepass changeit >/dev/null 2>&1 || true; \
      done; \
      break; \
    done

# The upstream builder scripts assume the sources sit in <src_dir>/quandenser.
COPY . /src/quandenser

# ProteoWizard and Boost. This is the slow step: it compiles Boost and the
# ProteoWizard data libraries with bjam and takes the better part of an hour.
# The upstream script fetches whatever ProteoWizard last built on the LabKey
# TeamCity server, so the version is not pinned; it is recorded in
# /opt/pwiz-version inside the image so a given image is at least traceable.
WORKDIR /build/tools
RUN /src/quandenser/ext/maracluster/admin/builders/install_proteowizard.sh /build/tools \
    && cp /build/tools/VERSION /opt/pwiz-version \
    && cat /opt/pwiz-version

# Rebuild the Boost include tree from one consistent version.
#
# install_proteowizard.sh leaves /build/tools/include/boost in a state that
# does not compile, for three compounding reasons.
#
# 1. It copies Boost's .ipp files from a hard-coded libraries/boost_1_76_0/
#    (its line 60), but the tarball it fetches from .lastSuccessful now ships
#    boost_1_86_0. That directory does not exist, rsync exits 23, and since the
#    script has no `set -e` the failure is discarded.
# 2. ProteoWizard bundles a *trimmed* Boost. Its 1.86 tree has no
#    boost/unordered/ at all, while the pinned Percolator includes
#    <boost/unordered/unordered_map.hpp> (ext/percolator/src/Scores.h:38).
# 3. The script finishes by unpacking boost_asio_1_18_2 from SourceForge and
#    rsyncing its whole boost/ directory over the bundled headers (lines
#    71-74), because Asio "is not included by the ProteoWizard boost tar but is
#    needed for maracluster". That held for the Boost 1.76 bundle it was
#    written against. It is now doubly wrong: the 1.86 bundle ships its own
#    Asio, and that SourceForge archive is a *standalone* Asio distribution
#    which carries its own copies of Boost.System, Boost.Config and others. So
#    the overlay does not add a missing library, it back-ports a
#    ten-release-old System and Asio over the current ones. The build fails
#    inside Boost's own headers, with thousands of lines of
#    "'BOOST_ASIO_CONSTEXPR' does not name a type" and "cannot convert
#    'boost::system::error_code' to 'int'".
#
# None of this is caused by the repair below: an untouched build of this
# repository fails in boost/asio/ as well.
#
# Patching the result file by file is a losing game, so discard the tree and
# lay it down once, in order: the bundled Boost as the authoritative base,
# then the matching upstream release to fill what ProteoWizard trimmed away,
# then ProteoWizard's own boost_aux additions. The obsolete Asio overlay is
# simply not reapplied.
RUN set -eux; \
    boost_dir="$(ls -d /build/tools/proteowizard/libraries/boost_1_*/ | head -n1)"; \
    boost_name="$(basename "${boost_dir%/}")"; \
    boost_dotted="$(echo "${boost_name#boost_}" | tr '_' '.')"; \
    inc=/build/tools/include/boost; \
    echo "Rebuilding the Boost include tree from ${boost_name}"; \
    rm -rf "${inc}"; \
    rsync -ap --include "*/" --include "*.hpp" --include "*.ipp" \
              --include "*.h" --exclude "*" \
          "${boost_dir}boost/" "${inc}/"; \
    wget --no-verbose -O /tmp/boost.tar.bz2 \
         "https://archives.boost.io/release/${boost_dotted}/source/${boost_name}.tar.bz2"; \
    mkdir -p /tmp/boostfull; \
    tar xf /tmp/boost.tar.bz2 -C /tmp/boostfull; \
    rsync -ap --ignore-existing \
              --include "*/" --include "*.hpp" --include "*.ipp" \
              --include "*.h" --exclude "*" \
          "/tmp/boostfull/${boost_name}/boost/" "${inc}/"; \
    rm -rf /tmp/boostfull /tmp/boost.tar.bz2; \
    if [ -d /build/tools/proteowizard/libraries/boost_aux/boost ]; then \
        rsync -ap --include "*/" --include "*.hpp" --include "*.ipp" \
                  --include "*.h" --exclude "*" \
              /build/tools/proteowizard/libraries/boost_aux/boost/ "${inc}/"; \
    fi; \
    test -f "${inc}/date_time/gregorian_calendar.ipp"; \
    test -f "${inc}/unordered/unordered_map.hpp"; \
    test -f "${inc}/asio/detail/reactive_socket_connect_op.hpp"; \
    grep -q "BOOST_ASIO_CONSTEXPR" "${inc}/asio/detail/config.hpp"; \
    echo "${boost_name}" > /opt/boost-version; \
    cat /opt/boost-version

# Quandenser itself, plus MaRaCluster and Percolator which build as
# subdirectories, plus the Dinosaur jar which CMake builds with Maven at
# configure time.
WORKDIR /build/quandenser
RUN cmake \
        -DTARGET_ARCH=amd64 \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr \
        -DBOOST_VARIANT_DO_NOT_USE_VARIADIC_TEMPLATES=ON \
        -DCMAKE_PREFIX_PATH=/build/tools \
        /src/quandenser \
    && make -j"${BUILD_JOBS}" \
    && make install

# ---------------------------------------------------------------------------
# Runtime
# ---------------------------------------------------------------------------
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.title="Quandenser" \
      org.opencontainers.image.description="Quantification by distillation for enhanced signals with error regulation" \
      org.opencontainers.image.source="https://github.com/statisticalbiotechnology/quandenser" \
      org.opencontainers.image.licenses="Apache-2.0"

COPY containers/quandenser/certs/ /tmp/extra-certs/
RUN set -eu; \
    if ls /tmp/extra-certs/*.crt >/dev/null 2>&1; then \
      mkdir -p /etc/ssl/certs; \
      cat /tmp/extra-certs/*.crt >> /etc/ssl/certs/ca-certificates.crt; \
      sed -i 's|http://|https://|g' /etc/apt/sources.list; \
    fi

# procps supplies `ps`, which Nextflow calls to collect per-task metrics.
# libgomp1 is the OpenMP runtime: Quandenser and MaRaCluster are compiled with
# -fopenmp, and the base image does not carry it (it arrives with gcc, which
# only exists in the builder stage).
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        libgomp1 \
        openjdk-8-jre-headless \
        procps \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/extra-certs

COPY --from=builder /usr/bin/quandenser /usr/bin/quandenser
COPY --from=builder /usr/share/java/ /usr/share/java/
COPY --from=builder /opt/pwiz-version /opt/pwiz-version
COPY --from=builder /opt/boost-version /opt/boost-version

# Smoke test. Deliberately not piped into anything: a pipeline would report the
# exit status of the last command and hide a binary that cannot even start,
# which is exactly how a missing libgomp.so.1 slips through.
#
# Note that Dinosaur's default Java heap is 24 GB (src/DinosaurIO.cpp:26) and
# that no environment variable can lower it, because the -Xmx is built into the
# command line Quandenser constructs. JAVA_TOOL_OPTIONS does not help. Callers
# must pass --dinosaur-memory; the Nextflow modules derive it from task.memory.
RUN set -eux; \
    quandenser --help 2>&1 | grep -q -i "usage\|quandenser"; \
    test -f /usr/share/java/Dinosaur-1.2.1.free.jar; \
    test -f /usr/share/java/advParams_dinosaur.txt; \
    test -f /usr/share/java/advParams_dinosaur_targeted.txt; \
    java -version

ENTRYPOINT []
CMD ["quandenser"]
