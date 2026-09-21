#!/bin/sh
#
# Pin the ProteoWizard source, and remove the floating downloads that
# ext/maracluster/admin/builders/install_proteowizard.sh would otherwise make.
#
# That script fetches whatever ProteoWizard's CI built most recently:
#
#   .../guestAuth/repository/download/bt81/.lastSuccessful/VERSION
#
# which is why this repository's build broke without a commit. Every hard-coded
# constant paired with that pointer (a Boost directory name, an Asio version, a
# C++ ABI flag) was correct when written and invalidated by a ProteoWizard
# release. Pinning removes the cause rather than the symptoms.
#
# The build is addressed the same way bioconda's own proteowizard recipe
# addresses it: by build id on the S3 bucket that TeamCity redirects to, with a
# checksum. Nothing here needs TeamCity to stay reachable.
#
# The bucket is not an archive. It keeps roughly the newest CI build and drops
# the rest, so a build id stops resolving within days: on 2026-09-21 both the
# id pinned here and the one current two hours earlier returned 403, leaving
# only that morning's build. Pinning therefore fixes *what* is built, not that
# it can still be fetched, and a from-scratch build needs a copy that does not
# depend on the bucket.
#
# Sources are tried in order: a file placed in containers/quandenser/cache/,
# then PWIZ_MIRROR_URL if set, then the bucket. Whichever answers, the checksum
# decides whether it is the pinned source.
#
# Usage: pin-proteowizard.sh <tools_dir> <build_id> <version> <sha256> <script> [cache_dir]
set -eu

tools_dir=$1
build_id=$2
version=$3
sha256=$4
script=$5
cache_dir=${6:-}

stem="pwiz-src-without-tv-$(echo "${version}" | tr ' ' '_')"
url="https://mc-tca-01.s3.us-west-2.amazonaws.com/ProteoWizard/bt81/${build_id}/${stem}.tar.bz2"

mkdir -p "${tools_dir}"
cd "${tools_dir}"

# install_proteowizard.sh derives the tarball name from this file, so writing
# it is what makes the rest of the script use our pinned copy.
printf '%s' "${version}" > VERSION

if [ -n "${cache_dir}" ] && [ -f "${cache_dir}/${stem}.tar.bz2" ]; then
  echo "Using the cached ProteoWizard ${version} from ${cache_dir}"
  cp "${cache_dir}/${stem}.tar.bz2" "${stem}.tar.bz2"
elif [ -n "${PWIZ_MIRROR_URL:-}" ]; then
  echo "Fetching pinned ProteoWizard ${version} from ${PWIZ_MIRROR_URL}"
  wget --no-verbose -O "${stem}.tar.bz2" "${PWIZ_MIRROR_URL}"
else
  echo "Fetching pinned ProteoWizard ${version} (build ${build_id})"
  if ! wget --no-verbose -O "${stem}.tar.bz2" "${url}"; then
    rm -f "${stem}.tar.bz2"
    echo "ERROR: ${url} could not be fetched." >&2
    echo "The bucket keeps only the newest build, so a pinned id stops resolving." >&2
    echo "Put ${stem}.tar.bz2 (sha256 ${sha256}) in containers/quandenser/cache/," >&2
    echo "or point PWIZ_MIRROR_URL at a copy, or re-pin all three PWIZ_ args." >&2
    exit 8
  fi
fi
echo "${sha256}  ${stem}.tar.bz2" | sha256sum -c -

# Neutralise the two fetches from the floating TeamCity pointer. The line that
# derives the filename from VERSION is deliberately left alone.
sed -i 's|^wget .*guestAuth.*$|: # source pinned by pin-proteowizard.sh|' "${script}"

# Neutralise the Boost.Asio overlay as well. It unpacks a standalone Asio 1.18.2
# distribution over ProteoWizard's bundled headers, which carries its own
# Boost.System and no longer matches the bundled Boost; the Dockerfile discards
# the result and rebuilds the header tree from one consistent version. Left in
# place it is not merely wasted work: the rsync is the script's last command, so
# a SourceForge outage fails the build over a file that is then thrown away.
sed -i \
    -e 's|^wget .*sourceforge.*$|: # asio overlay removed by pin-proteowizard.sh|' \
    -e 's|^tar -xzf boost_asio.tar.gz$|: # asio overlay removed by pin-proteowizard.sh|' \
    -e 's|^rsync .*boost_asio_1_18_2.*$|: # asio overlay removed by pin-proteowizard.sh|' \
    "${script}"

# Fail loudly if upstream reshapes these lines, rather than silently building
# against a floating source again.
for pattern in guestAuth sourceforge boost_asio_1_18_2; do
    if grep -q "${pattern}" "${script}"; then
        echo "ERROR: '${pattern}' still present in ${script}; the pin did not apply" >&2
        grep -n "${pattern}" "${script}" >&2
        exit 1
    fi
done

echo "ProteoWizard pinned to ${version} (build ${build_id}); floating downloads removed"
