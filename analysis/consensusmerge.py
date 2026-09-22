#!/usr/bin/env python3
"""
Python ports of MaRaCluster's two consensus merges.

`merge_thomson` is MSClusterMerge: peaks are joined within a fixed number of
Thomson, first inside each spectrum (0.34 Th) and then across the cluster
(0.22 Th). `merge_ppm` is PpmConsensusMerge: peaks are grouped on a ppm scale,
no two peaks of one spectrum ever join, and the peak cap can be switched off.

They exist here because the analysis scripts build their own merged spectra and
have to be able to build them the way the pipeline now can, without a C++
build. `consensusmerge_check.py` compiles the two C++ sources standalone and
checks that these functions reproduce them peak for peak.

Both take a cluster as [(mz, intensity), ...], one entry per member spectrum,
and return (mz, intensity) sorted by m/z. Members are scaled to a total
intensity of 1000 first, as MSFileMerger does, so the caller's normalization
does not change the result.
"""
import numpy as np

# PpmConsensusMerge defaults
PPM_SIGMA = 2.0       # fragment m/z scatter between two spectra, in ppm
NUM_SIGMAS = 4.0      # peaks further apart than this many sigmas are separate
MAX_PEAKS = 160       # most intense peaks to keep; 0 keeps all of them

# MSClusterMerge constants, inherited from MS-Cluster
MASS_TO_INT_RATIO = np.float32(10000.0)
FRAGMENT_TOLERANCE = np.float32(0.34)
ISO_TOLERANCE = np.float32(np.float32(0.1) +
                           (FRAGMENT_TOLERANCE - np.float32(0.1)) *
                           np.float32(0.5))
THOMSON_MAX_PEAKS = 160   # not configurable in MSClusterMerge


class PeakWeightTable:
    """Binomial occurrence weight, as in MSClusterMerge's table.

    A group seen in k of n members is discounted by P(Binomial(n, p) < k), the
    chance that noise alone would have produced fewer occurrences. Computed in
    float32 because the C++ table is a vector<float> and its 0.99 cutoff is
    taken on that value.
    """

    def __init__(self, max_n=64, p=float(np.float32(0.15))):
        # p is a float in the C++ signature, so 0.15 is not exactly 0.15 there
        # either; taking the same rounded value keeps the tables identical.
        self.max_n = max_n
        self.weights = [[np.float32(0.0)] for _ in range(max_n + 1)]
        self.k_with_value_one = [0] * (max_n + 1)
        for n in range(1, max_n + 1):
            w = [np.float32(0.0)]
            cdf = np.float32(0.0)
            coeff = 1.0
            power = (1.0 - p) ** n
            ratio = p / (1.0 - p)
            first_one = n + 1
            for k in range(0, n + 1):
                if k == 0:
                    cdf = np.float32(power)
                else:
                    coeff *= (n - k + 1.0) / k
                    power *= ratio
                    cdf = np.float32(cdf + np.float32(coeff * power))
                if cdf > np.float32(0.99):
                    first_one = k
                    break
                w.append(cdf)
            self.weights[n] = w
            self.k_with_value_one[n] = first_one

    def get_weight(self, k, n):
        while n > self.max_n:
            n >>= 1
            k >>= 1
        if k >= self.k_with_value_one[n]:
            return 1.0
        return float(self.weights[n][k])


_TABLE = PeakWeightTable()


def _normalized(cluster):
    """Each member's positive peaks, scaled to a total intensity of 1000."""
    out = []
    for mz, intensity in cluster:
        mz = np.asarray(mz, np.float64)
        intensity = np.asarray(intensity, np.float64)
        keep = (intensity > 0.0) & (mz > 0.0)
        total = intensity[keep].sum()
        if total <= 0.0:
            out.append((np.empty(0), np.empty(0)))
            continue
        # the normalization factor is a float in SpectrumHandler, and the
        # Thomson merge keeps the difference all the way to its output
        scale = float(np.float32(1000.0 / total))
        o = np.argsort(mz[keep], kind='stable')
        out.append((mz[keep][o], intensity[keep][o] * scale))
    return out


# --- PpmConsensusMerge ----------------------------------------------------

def _pool_peaks(cluster):
    log_mz, intensity, member = [], [], []
    for i, (mz, it) in enumerate(cluster):
        if mz.size == 0:
            continue
        total = it.sum()
        if total <= 0.0:
            continue
        # poolPeaks scales each member to 1000 again, in double, on top of
        # whatever normalization MSFileMerger already applied
        log_mz.append(np.log(mz) * 1e6)
        intensity.append(it * (1000.0 / total))
        member.append(np.full(mz.size, i, np.int64))
    if not log_mz:
        return np.empty(0), np.empty(0), np.empty(0, np.int64)
    log_mz = np.concatenate(log_mz)
    intensity = np.concatenate(intensity)
    member = np.concatenate(member)
    o = np.argsort(log_mz, kind='stable')
    return log_mz[o], intensity[o], member[o]


def _segment_component(log_mz, member, begin, end, num_members, next_group,
                       group_idx, ppm_sigma, num_sigmas):
    """Split one component so that no group holds two peaks of one spectrum.

    Dynamic programming over m/z order: within-group scatter in units of the
    single-peak variance, plus a constant num_sigmas**2 per group.
    """
    n = end - begin
    seen = np.full(num_members, -1, np.int64)
    duplicate = False
    for i in range(begin, end):
        if seen[member[i]] == 0:
            duplicate = True
            break
        seen[member[i]] = 0

    if not duplicate:
        group_idx[begin:end] = next_group
        return next_group + 1

    gap = num_sigmas * ppm_sigma
    peak_variance = 0.5 * ppm_sigma * ppm_sigma
    penalty = num_sigmas * num_sigmas

    best_cost = np.full(n + 1, np.inf)
    best_start = np.zeros(n + 1, np.int64)
    best_cost[0] = 0.0
    seen[:] = -1

    for i in range(1, n + 1):
        total = sq = 0.0
        count = 0
        for j in range(i - 1, -1, -1):
            m = member[begin + j]
            # a group holds at most one peak per spectrum, and spans at most
            # twice the gap that would have cut it in the first place
            if seen[m] == i:
                break
            if log_mz[begin + i - 1] - log_mz[begin + j] > 2.0 * gap:
                break
            seen[m] = i
            total += log_mz[begin + j]
            sq += log_mz[begin + j] * log_mz[begin + j]
            count += 1
            scatter = sq - total * total / count
            cost = scatter / peak_variance + penalty + best_cost[j]
            if cost < best_cost[i]:
                best_cost[i] = cost
                best_start[i] = j

    starts = []
    i = n
    while i > 0:
        starts.append(int(best_start[i]))
        i = int(best_start[i])
    for s in range(len(starts) - 1, -1, -1):
        segment_end = n if s == 0 else starts[s - 1]
        group_idx[begin + starts[s]:begin + segment_end] = next_group
        next_group += 1
    return next_group


def _group_peaks(log_mz, member, num_members, ppm_sigma, num_sigmas):
    group_idx = np.zeros(log_mz.size, np.int64)
    if log_mz.size == 0:
        return group_idx, 0
    gap = num_sigmas * ppm_sigma
    cuts = np.flatnonzero(np.diff(log_mz) > gap) + 1
    next_group = 0
    begin = 0
    for end in list(cuts) + [log_mz.size]:
        next_group = _segment_component(log_mz, member, int(begin), int(end),
                                        num_members, next_group, group_idx,
                                        ppm_sigma, num_sigmas)
        begin = end
    return group_idx, next_group


def merge_ppm(cluster, ppm_sigma=PPM_SIGMA, num_sigmas=NUM_SIGMAS,
              max_peaks=MAX_PEAKS):
    """PpmConsensusMerge: grouping on a ppm scale, one peak per member."""
    num_members = len(cluster)
    if num_members == 0:
        return np.empty(0), np.empty(0)

    log_mz, intensity, member = _pool_peaks(_normalized(cluster))
    if log_mz.size == 0:
        return np.empty(0), np.empty(0)

    group_idx, num_groups = _group_peaks(log_mz, member, num_members,
                                         ppm_sigma, num_sigmas)
    if num_groups == 0:
        return np.empty(0), np.empty(0)

    sum_intensity = np.bincount(group_idx, weights=intensity,
                                minlength=num_groups)
    sum_weighted = np.bincount(group_idx, weights=intensity * log_mz,
                               minlength=num_groups)
    # one peak per member per group, so this counts members
    counts = np.bincount(group_idx, minlength=num_groups)

    keep = sum_intensity > 0.0
    mz = np.exp(sum_weighted[keep] / sum_intensity[keep] * 1e-6)
    weights = np.array([_TABLE.get_weight(int(k), num_members)
                        for k in counts[keep]])
    it = sum_intensity[keep] * weights

    if max_peaks > 0 and mz.size > max_peaks:
        sel = np.argpartition(it, -max_peaks)[-max_peaks:]
        mz, it = mz[sel], it[sel]
    o = np.argsort(mz, kind='stable')
    return mz[o], it[o]


# --- MSClusterMerge -------------------------------------------------------

def _mass_to_int(mass):
    return int(np.float32(MASS_TO_INT_RATIO * np.float32(mass)))


def _bin_member(mz, intensity):
    """Join peaks of one spectrum within FRAGMENT_TOLERANCE, in m/z order."""
    if mz.size == 0:
        return [], [], []
    max_proximity = _mass_to_int(FRAGMENT_TOLERANCE)
    bins = [_mass_to_int(mz[0])]
    mzs = [float(mz[0])]
    its = [float(intensity[0])]
    for i in range(1, mz.size):
        cur = _mass_to_int(mz[i])
        if cur - bins[-1] < max_proximity:
            total = its[-1] + float(intensity[i])
            ratio = its[-1] / total
            new_mass = ratio * mzs[-1] + (1.0 - ratio) * float(mz[i])
            mzs[-1] = new_mass
            bins[-1] = _mass_to_int(new_mass)
            its[-1] = total
        else:
            bins.append(cur)
            mzs.append(float(mz[i]))
            its.append(float(intensity[i]))
    return bins, mzs, its


def merge_thomson(cluster):
    """MSClusterMerge: the merge the pipeline uses unless told otherwise."""
    num_members = len(cluster)
    if num_members == 0:
        return np.empty(0), np.empty(0)

    pooled = []
    for i, (mz, it) in enumerate(_normalized(cluster)):
        bins, mzs, its = _bin_member(mz, it)
        pooled += [(b, m, v, i) for b, m, v in zip(bins, mzs, its) if v > 0.0]
    if not pooled:
        return np.empty(0), np.empty(0)
    pooled.sort(key=lambda p: (p[1], p[2]))

    max_proximity = _mass_to_int(ISO_TOLERANCE)
    bins = [pooled[0][0]]
    mzs = [pooled[0][1]]
    its = [pooled[0][2]]
    members = {pooled[0][3]}
    counts = []
    for b, m, v, member in pooled[1:]:
        if b - bins[-1] < max_proximity:
            total = its[-1] + v
            ratio = its[-1] / total
            new_mass = ratio * mzs[-1] + (1.0 - ratio) * m
            mzs[-1] = new_mass
            bins[-1] = _mass_to_int(new_mass)
            its[-1] = total
        else:
            counts.append(len(members))
            members = set()
            bins.append(b)
            mzs.append(m)
            its.append(v)
        members.add(member)
    counts.append(len(members))

    # discount peaks that have only a few copies
    mz = np.array(mzs)
    it = np.array([v * _TABLE.get_weight(c, num_members)
                   for v, c in zip(its, counts)])

    if mz.size > THOMSON_MAX_PEAKS:
        sel = np.argpartition(it, -THOMSON_MAX_PEAKS)[-THOMSON_MAX_PEAKS:]
        mz, it = mz[sel], it[sel]
    o = np.argsort(mz, kind='stable')
    return mz[o], it[o]
