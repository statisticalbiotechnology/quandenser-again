#!/usr/bin/env python3
"""
Check `consensusmerge.py` against the C++ it is a port of.

Compiles `PpmConsensusMerge.cpp` and `MSClusterMerge.cpp` out of
`ext/maracluster/src` on their own, against stub headers standing in for the
ProteoWizard-dependent ones, runs the ppm merge's unit test, and compares both
merges with the Python ones on random clusters built to hit the cases they
differ on: replicate peaks a few ppm apart, isotope pairs a few tenths of a
Thomson apart, two fragments inside one gap, and peaks seen in one member only.

Only the normalization to a total intensity of 1000, five lines of
SpectrumHandler, is transcribed into the harness; the merges themselves are the
repository's own sources, compiled unmodified. Needs g++ and the Boost headers
(apt-get install libboost-dev), which MSClusterMerge.cpp includes.

    python3 consensusmerge_check.py [trials]

Reports the largest disagreement in m/z, in ppm, and exits non-zero on any.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

import consensusmerge

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   '..', 'ext', 'maracluster', 'src')

STUBS = {
    'MZIntensityPair.h': """
#ifndef STUB_MZINTENSITYPAIR_H_
#define STUB_MZINTENSITYPAIR_H_
/* the real header reaches these through pwiz's MSData.hpp */
#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <vector>
namespace maracluster {
struct MZIntensityPair {
  double mz, intensity, multiplicity;
  MZIntensityPair() : mz(0.0), intensity(0.0), multiplicity(1.0) {}
  MZIntensityPair(double mz_, double intensity_)
      : mz(mz_), intensity(intensity_), multiplicity(1.0) {}
};
}
#endif
""",
    'MassChargeCandidate.h': """
#ifndef STUB_MASSCHARGECANDIDATE_H_
#define STUB_MASSCHARGECANDIDATE_H_
namespace maracluster {
struct MassChargeCandidate {
  MassChargeCandidate(unsigned int charge_, double precMz_, double mass_)
      : charge(charge_), precMz(precMz_), mass(mass_) {}
  unsigned int charge;
  double precMz, mass;
  inline static bool lessChargeMass(const MassChargeCandidate& a,
                                    const MassChargeCandidate& b) {
    return (a.charge < b.charge) || ((a.charge == b.charge) && (a.mass < b.mass));
  }
};
}
#endif
""",
    'SpectrumHandler.h': """
#ifndef STUB_SPECTRUMHANDLER_H_
#define STUB_SPECTRUMHANDLER_H_
#include "MZIntensityPair.h"
#include "MassChargeCandidate.h"
namespace maracluster {
class SpectrumHandler {
 public:
  static const double PROTON_MASS;
  inline static bool lessMZ(const MZIntensityPair& a, const MZIntensityPair& b) {
    return (a.mz < b.mz) || (a.mz == b.mz && a.intensity < b.intensity);
  }
  inline static bool greaterIntensity(const MZIntensityPair& a, const MZIntensityPair& b) {
    return (a.intensity > b.intensity);
  }
  inline static double calcPrecMz(double mass, unsigned int charge) {
    return (mass + PROTON_MASS * (charge - 1)) / charge;
  }
};
}
#endif
""",
    'BinSpectra.h': """
#ifndef STUB_BINSPECTRA_H_
#define STUB_BINSPECTRA_H_
#include <vector>
#include <iostream>
#include <cstring>
#include "MZIntensityPair.h"
#include "SpectrumHandler.h"
namespace maracluster {
class BinSpectra {
 public:
  static inline double getMZ(unsigned int bin) { return static_cast<double>(bin); }
};
struct BinnedMZIntensityPair : public MZIntensityPair {
 public:
  unsigned int binIdx;
  BinnedMZIntensityPair(unsigned int binIdx_, double mz_, double intensity_,
                        double multiplicity_) {
    binIdx = binIdx_; mz = mz_; intensity = intensity_; multiplicity = multiplicity_;
  }
  BinnedMZIntensityPair(unsigned int binIdx_, double mz_, double intensity_) {
    binIdx = binIdx_; mz = mz_; intensity = intensity_; multiplicity = 1.0;
  }
  BinnedMZIntensityPair() { binIdx = 0; mz = 0.0; intensity = 0.0; multiplicity = 1.0; }
};
}
#endif
""",
    'main.cpp': """
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include "PpmConsensusMerge.h"
#include "MSClusterMerge.h"
using namespace maracluster;
const double SpectrumHandler::PROTON_MASS = 1.00727646677;

static void readCluster(std::vector<std::vector<MZIntensityPair> >& cluster) {
  size_t numMembers = 0;
  if (std::scanf("%zu", &numMembers) != 1) exit(2);
  cluster.resize(numMembers);
  for (size_t i = 0; i < numMembers; ++i) {
    size_t n = 0;
    if (std::scanf("%zu", &n) != 1) exit(2);
    for (size_t j = 0; j < n; ++j) {
      double mz, intensity;
      if (std::scanf("%lf %lf", &mz, &intensity) != 2) exit(2);
      cluster[i].push_back(MZIntensityPair(mz, intensity));
    }
  }
}

static void writeSpectrum(const std::vector<MZIntensityPair>& merged) {
  std::printf("%zu\\n", merged.size());
  for (size_t i = 0; i < merged.size(); ++i) {
    std::printf("%.17g %.17g\\n", merged[i].mz, merged[i].intensity);
  }
}

/* the five lines of SpectrumHandler::normalizeIntensitiesMSCluster that
   MSFileMerger applies before either merge */
static void normalize(std::vector<MZIntensityPair>& mziPairs) {
  double total = 0.0;
  for (size_t i = 0; i < mziPairs.size(); ++i) total += mziPairs[i].intensity;
  if (total <= 0.0) return;
  const float scale = 1000.0 / total;
  for (size_t i = 0; i < mziPairs.size(); ++i) mziPairs[i].intensity *= scale;
}

int main(int argc, char** argv) {
  const char* mode = (argc > 1) ? argv[1] : "ppm";
  if (std::strcmp(mode, "unittest") == 0) {
    bool ok = PpmConsensusMerge::mergeUnitTest();
    std::cout << (ok ? "unit test PASSED" : "unit test FAILED") << std::endl;
    return ok ? 0 : 1;
  }
  std::vector<std::vector<MZIntensityPair> > cluster;
  std::vector<MZIntensityPair> merged;
  if (std::strcmp(mode, "thomson") == 0) {
    MSClusterMerge::init();
    readCluster(cluster);
    std::vector<std::vector<BinnedMZIntensityPair> > binned(cluster.size());
    for (size_t i = 0; i < cluster.size(); ++i) {
      normalize(cluster[i]);
      MSClusterMerge::binMZIntensityPairs(cluster[i], binned[i]);
    }
    std::vector<BinnedMZIntensityPair> mergedBinned;
    MSClusterMerge::merge(binned, mergedBinned);
    MSClusterMerge::unbinMZIntensityPairs(mergedBinned, merged);
    writeSpectrum(merged);
    return 0;
  }
  if (argc > 4) {
    PpmConsensusMerge::ppmSigma_ = static_cast<float>(atof(argv[2]));
    PpmConsensusMerge::numSigmas_ = static_cast<float>(atof(argv[3]));
    PpmConsensusMerge::maxPeaks_ = static_cast<unsigned int>(atoi(argv[4]));
  }
  PpmConsensusMerge::init();
  readCluster(cluster);
  for (size_t i = 0; i < cluster.size(); ++i) normalize(cluster[i]);
  PpmConsensusMerge::merge(cluster, merged);
  writeSpectrum(merged);
  return 0;
}
""",
}


def build(workdir):
    """Compile both merges against the stubs. Returns the binary path."""
    for name, text in STUBS.items():
        with open(os.path.join(workdir, name), 'w') as fh:
            fh.write(text)
    for name in ('PpmConsensusMerge.h', 'PpmConsensusMerge.cpp',
                 'MSClusterMerge.h', 'MSClusterMerge.cpp'):
        with open(os.path.join(SRC, name)) as src, \
                open(os.path.join(workdir, name), 'w') as dst:
            dst.write(src.read())
    binary = os.path.join(workdir, 'consensusmerge')
    subprocess.run(['g++', '-O2', '-o', binary, 'main.cpp',
                    'PpmConsensusMerge.cpp', 'MSClusterMerge.cpp'],
                   cwd=workdir, check=True)
    return binary


def run_cpp(binary, cluster, args):
    lines = ['%d' % len(cluster)]
    for mz, it in cluster:
        lines.append('%d' % len(mz))
        lines += ['%.17g %.17g' % (a, b) for a, b in zip(mz, it)]
    out = subprocess.run([binary] + [str(a) for a in args],
                         input='\n'.join(lines) + '\n', check=True,
                         capture_output=True, text=True).stdout.split('\n')
    n = int(out[0])
    mz = np.array([float(out[1 + i].split()[0]) for i in range(n)])
    it = np.array([float(out[1 + i].split()[1]) for i in range(n)])
    return mz, it


def random_cluster(rng):
    """A cluster with the peak configurations the merges have to get right."""
    num_members = int(rng.integers(2, 9))
    members = [[] for _ in range(num_members)]

    for _ in range(int(rng.integers(3, 12))):       # replicate fragment
        mz0 = float(rng.uniform(300.0, 1800.0))
        base = float(rng.uniform(10.0, 1000.0))
        for i in range(num_members):
            if rng.random() < 0.25:                 # missing in some members
                continue
            members[i].append((mz0 * (1.0 + rng.normal(0.0, 2e-6)),
                               base * float(rng.uniform(0.5, 2.0))))

    for _ in range(int(rng.integers(0, 4))):        # isotope pair, one charge
        mz0 = float(rng.uniform(300.0, 1800.0))
        spacing = 1.00286 / float(rng.integers(1, 5))
        base = float(rng.uniform(10.0, 1000.0))
        for i in range(num_members):
            members[i].append((mz0 * (1.0 + rng.normal(0.0, 2e-6)), base))
            members[i].append(((mz0 + spacing) * (1.0 + rng.normal(0.0, 2e-6)),
                               base * 0.5))

    for _ in range(int(rng.integers(0, 4))):        # two fragments in one gap
        mz0 = float(rng.uniform(300.0, 1800.0))
        delta = mz0 * float(rng.uniform(4e-6, 8e-6))
        for i in range(num_members):
            if rng.random() < 0.3:
                continue
            members[i].append((mz0 * (1.0 + rng.normal(0.0, 2e-6)),
                               float(rng.uniform(10.0, 1000.0))))
            members[i].append(((mz0 + delta) * (1.0 + rng.normal(0.0, 2e-6)),
                               float(rng.uniform(10.0, 1000.0))))

    for _ in range(int(rng.integers(0, 15))):       # noise, one member only
        i = int(rng.integers(0, num_members))
        members[i].append((float(rng.uniform(300.0, 1800.0)),
                           float(rng.uniform(1.0, 100.0))))

    cluster = []
    for peaks in members:
        peaks.sort()
        cluster.append((np.array([p[0] for p in peaks]),
                        np.array([p[1] for p in peaks])))
    return cluster


def compare(name, py, cpp, trial, worst):
    py_mz, py_it = py
    c_mz, c_it = cpp
    if py_mz.size != c_mz.size:
        print('trial %d %s: %d peaks in python, %d in C++'
              % (trial, name, py_mz.size, c_mz.size))
        return worst, 1
    if py_mz.size == 0:
        return worst, 0
    d_mz = np.abs(py_mz - c_mz) / c_mz * 1e6
    d_it = np.abs(py_it - c_it) / np.maximum(c_it, 1e-12)
    worst = (max(worst[0], float(d_mz.max())), max(worst[1], float(d_it.max())))
    if d_mz.max() > 1e-6 or d_it.max() > 1e-9:
        print('trial %d %s: max %.3g ppm in m/z, %.3g relative in intensity'
              % (trial, name, d_mz.max(), d_it.max()))
        return worst, 1
    return worst, 0


def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    rng = np.random.default_rng(20260922)
    workdir = tempfile.mkdtemp(prefix='consensusmerge-')
    binary = build(workdir)
    print(subprocess.run([binary, 'unittest'], check=True, capture_output=True,
                         text=True).stdout.strip())

    worst_ppm = (0.0, 0.0)
    worst_thomson = (0.0, 0.0)
    failures = 0
    for trial in range(trials):
        ppm_sigma = float(rng.choice([1.0, 2.0, 5.0, 20.0]))
        num_sigmas = float(rng.choice([2.0, 4.0]))
        max_peaks = int(rng.choice([0, 10, 160]))
        cluster = random_cluster(rng)

        worst_ppm, bad = compare(
            'ppm',
            consensusmerge.merge_ppm(cluster, ppm_sigma, num_sigmas, max_peaks),
            run_cpp(binary, cluster, ['ppm', ppm_sigma, num_sigmas, max_peaks]),
            trial, worst_ppm)
        failures += bad

        worst_thomson, bad = compare(
            'thomson',
            consensusmerge.merge_thomson(cluster),
            run_cpp(binary, cluster, ['thomson']),
            trial, worst_thomson)
        failures += bad

    print('%d trials' % trials)
    print('  ppm      largest disagreement %.3g ppm in m/z, %.3g relative in '
          'intensity' % worst_ppm)
    print('  thomson  largest disagreement %.3g ppm in m/z, %.3g relative in '
          'intensity' % worst_thomson)
    print('FAILED on %d comparisons' % failures if failures else 'identical')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
