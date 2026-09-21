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
#
# The source is pinned to one build. Upstream fetches whatever ProteoWizard's
# CI built most recently, which is the single reason this repository's build
# broke with no commit to it: every hard-coded constant paired with that
# floating pointer went stale on someone else's release schedule. Change these
# three values together, deliberately, and expect to revisit the Boost repair
# below when you do.
#
# To move to a newer ProteoWizard, take the build id from the S3 redirect of
#   https://teamcity.labkey.org/guestAuth/repository/download/bt81/.lastSuccessful/VERSION
# and the checksum from the tarball it names.
ARG PWIZ_BUILD_ID=4176587
ARG PWIZ_VERSION="3 0 26258 7993a4e"
ARG PWIZ_SHA256=d36afbb6736cc44ec5bbb38a482de9f2d4c3d962c2d6b059f423bccbff63e10d

COPY containers/quandenser/pin-proteowizard.sh /usr/local/bin/pin-proteowizard.sh

# A copy of the pinned tarball, if one was put there. ProteoWizard's bucket
# keeps only the newest CI build, so the pinned build id stops resolving after
# a few days and a from-scratch build then fails at this step unless it finds
# the source here (or at PWIZ_MIRROR_URL). The directory normally holds
# nothing but .gitkeep.
COPY containers/quandenser/cache/ /tmp/pwiz-cache/
ARG PWIZ_MIRROR_URL=""

WORKDIR /build/tools
RUN PWIZ_MIRROR_URL="${PWIZ_MIRROR_URL}" /usr/local/bin/pin-proteowizard.sh \
        /build/tools \
        "${PWIZ_BUILD_ID}" \
        "${PWIZ_VERSION}" \
        "${PWIZ_SHA256}" \
        /src/quandenser/ext/maracluster/admin/builders/install_proteowizard.sh \
        /tmp/pwiz-cache \
    && /src/quandenser/ext/maracluster/admin/builders/install_proteowizard.sh /build/tools \
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
# 3. The script used to finish by unpacking boost_asio_1_18_2 from SourceForge
#    over the bundled headers, because Asio "is not included by the
#    ProteoWizard boost tar but is needed for maracluster". That held for the
#    Boost 1.76 bundle it was written against. It became doubly wrong: the 1.86
#    bundle ships its own Asio, and that SourceForge archive is a *standalone*
#    Asio distribution carrying its own Boost.System and Boost.Config, so the
#    overlay back-ported a ten-release-old System next to the current Asio and
#    the build failed inside Boost's own headers. pin-proteowizard.sh now
#    removes that overlay before the script runs.
#
# None of this is caused by the repair below: an untouched build of this
# repository fails in boost/asio/ as well.
#
# Patching the result file by file is a losing game, so discard the tree and
# lay it down once, in order: the bundled Boost as the authoritative base,
# then the matching upstream release to fill what ProteoWizard trimmed away,
# then ProteoWizard's own boost_aux additions.
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
    test -f "${inc}/version.hpp"; \
    echo "${boost_name}" > /opt/boost-version; \
    cat /opt/boost-version

# Give zstd the same unsuffixed symlink the Boost libraries get.
#
# ProteoWizard's bjam names its output libzstd-gcc11-mt-s.a, and
# install_proteowizard.sh creates plain-name symlinks only for the Boost
# components (its lines 66-68), because when it was written nothing else
# needed one. Current ProteoWizard compresses binary data with zstd, so
# libpwiz_data_msdata_core now references ZSTD_compress and friends, and with
# no libzstd.a there is nothing for -lzstd to find.
RUN set -eux; \
    zstd_lib="$(ls /build/tools/lib/libzstd-*.a | head -n1)"; \
    ln -sf "${zstd_lib}" /build/tools/lib/libzstd.a; \
    ls -la /build/tools/lib/libzstd.a

# Quandenser itself, plus MaRaCluster and Percolator which build as
# subdirectories, plus the Dinosaur jar which CMake builds with Maven at
# configure time.
# Align the C++ string ABI with the one ProteoWizard is actually built with.
#
# ext/maracluster/CommonCMake.txt sets _GLIBCXX_USE_CXX11_ABI=0
# unconditionally (its lines 106-110), which selects the pre-C++11
# std::basic_string. That matched the ProteoWizard of the day. bjam now builds
# ProteoWizard with the current default, so its archives define
# std::__cxx11::basic_string while everything compiled here references the old
# one. The mangled names differ, so the link fails with undefined references
# to functions the archives demonstrably contain: a symbol dump shows
# pwiz::msdata::MSDataFile::MSDataFile defined, taking
# std::__cxx11::basic_string, against an undefined reference taking
# std::string.
#
# Everything else in this image is compiled from source in one pass, so moving
# to the current ABI makes the whole build self-consistent, whereas keeping
# ABI=0 would mean rebuilding ProteoWizard and Boost to match it.
RUN sed -i 's/_GLIBCXX_USE_CXX11_ABI=0/_GLIBCXX_USE_CXX11_ABI=1/g' \
        /src/quandenser/ext/maracluster/CommonCMake.txt \
    && grep -n "_GLIBCXX_USE_CXX11_ABI" /src/quandenser/ext/maracluster/CommonCMake.txt

# -include boost/assert.hpp compensates for a transitive include that Boost
# withdrew. MaRaCluster uses BOOST_ASSERT in PvalueCalculator.cpp without
# including <boost/assert.hpp>; it used to arrive by way of some other Boost
# header, and with 1.86 it no longer does, so the file stops compiling. This is
# the one place where the age gap between the pinned submodules and the Boost
# that current ProteoWizard bundles shows up in code this project owns, rather
# than in the header tree. Forcing the include is a smaller change than
# carrying a patch against a submodule, and boost/assert.hpp is small and
# idempotent, so applying it to every translation unit costs nothing.
#
# The path has to be absolute. CMake's own compiler probes compile with
# CMAKE_CXX_FLAGS but without the project's include directories, so a bare
# "-include boost/assert.hpp" fails during configure, long before anything in
# this project is built.
#
# CMAKE_CXX_STANDARD_LIBRARIES appends to the end of the link line, which is
# where the remaining two problems are fixed. These static archives reference
# each other, and ld resolves a link line left to right in one pass, so
# whichever order FindProteoWizard happens to produce leaves symbols such as
# pwiz::msdata::MSDataFile::MSDataFile undefined even though the archive
# defining them is right there. Repeating the set inside --start-group makes
# ld iterate until nothing more resolves. The group also pulls in zstd, which
# nothing in the CMake files knows about because ProteoWizard only started
# compressing binary data with it after these build files were written.
WORKDIR /build/quandenser
RUN cmake \
        -DTARGET_ARCH=amd64 \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr \
        -DBOOST_VARIANT_DO_NOT_USE_VARIADIC_TEMPLATES=ON \
        -DCMAKE_PREFIX_PATH=/build/tools \
        -DCMAKE_CXX_FLAGS="-include /build/tools/include/boost/assert.hpp" \
        -DCMAKE_CXX_STANDARD_LIBRARIES="-Wl,--start-group \
            -L/build/tools/lib \
            -lpwiz_data_msdata -lpwiz_data_msdata_core -lpwiz_data_common \
            -lpwiz_utility_misc -lpwiz_utility_minimxml -lSHA1 \
            -lboost_serialization -lboost_iostreams -lboost_filesystem \
            -lboost_system -lboost_thread -lboost_chrono -lboost_regex \
            -lzstd -lzlib \
            -Wl,--end-group" \
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
# Only the three files Quandenser looks for. /usr/share/java in the builder
# also holds every jar Maven brought in, which is some 150 files and tens of
# megabytes that this image has no use for.
COPY --from=builder /usr/share/java/Dinosaur-1.2.1.free.jar /usr/share/java/
COPY --from=builder /usr/share/java/advParams_dinosaur.txt /usr/share/java/
COPY --from=builder /usr/share/java/advParams_dinosaur_targeted.txt /usr/share/java/
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
