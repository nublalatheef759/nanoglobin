"""
SV identification against the IthaCNVs catalogue (databases/cnvs.csv).

Matching requires BOTH:

  1. Reciprocal overlap >= MIN_RECIPROCAL in both directions. One-directional
     containment is not sufficient: a 92bp deletion sits entirely inside the
     Hb Lepore window but is not Hb Lepore, and a true -a3.7 sits entirely
     inside the --SEA window but is not --SEA.

  2. Size concordance within SIZE_TOL. This is what separates -a3.7 (~3.8kb)
     from -a4.2 (~4.3kb), whose windows overlap substantially.

Breakpoint uncertainty is handled explicitly. IthaCNVs reports 67/258 entries as
ranges, because those junctions were characterised by gap-PCR or MLPA and never
sequenced. For those, the catalogue gives an interval the breakpoint lies within.
We take the most favourable position consistent with that stated interval, rather
than collapsing it to a point and inventing precision the literature lacks.
"""

import csv
import os

MIN_RECIPROCAL = 0.5
SIZE_TOL = 0.10

_DEFAULT_DB = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "databases", "cnvs.csv",
)

# A candidate DEL should be allowed to match a catalogued DELINS: those are
# deletions with a few bases inserted at the junction, which a caller reports
# as a plain deletion. COMPLEX and untyped entries are not matchable.
_COMPATIBLE = {
    "DEL": {"DEL", "DELINS"},
    "DUP": {"DUP"},
}

_CATALOGUE = None


def load_catalogue(path=None):
    """Load databases/cnvs.csv. Returns a list of dicts; [] if the file is absent."""
    path = path or _DEFAULT_DB
    if not os.path.exists(path):
        return []
    out = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            try:
                entry = {
                    "name": row["name"],
                    "chrom": row["chrom"],
                    "start_min": int(row["start_min"]),
                    "start_max": int(row["start_max"]),
                    "end_min": int(row["end_min"]),
                    "end_max": int(row["end_max"]),
                    "type": row["type"],
                    "ithaID": row.get("ithaID", ""),
                    "exact": row.get("exact", ""),
                }
            except (KeyError, ValueError, TypeError):
                continue
            if not entry["chrom"] or entry["type"] not in ("DEL", "DELINS", "DUP"):
                continue
            out.append(entry)
    return out


def _catalogue():
    global _CATALOGUE
    if _CATALOGUE is None:
        _CATALOGUE = load_catalogue()
    return _CATALOGUE


def _clamp(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)


def identify_sv(chrom, pos, svtype, svlen):
    """Return a catalogued variant name, or unknown_<TYPE>_<LEN>bp."""
    pos = int(pos)
    svlen = abs(int(svlen)) if svlen else 0
    if svlen == 0:
        return "unknown_%s_0bp" % svtype
    cand_start, cand_end = pos, pos + svlen

    compatible = _COMPATIBLE.get(svtype)
    if compatible is None:
        return "unknown_%s_%dbp" % (svtype, svlen)

    best_name, best_score, best_itha = None, 0.0, ""

    for e in _catalogue():
        if e["chrom"] != chrom or e["type"] not in compatible:
            continue

        # Most favourable breakpoints consistent with the catalogued interval.
        k_start = _clamp(cand_start, e["start_min"], e["start_max"])
        k_end = _clamp(cand_end, e["end_min"], e["end_max"])
        known_len = k_end - k_start
        if known_len <= 0:
            continue

        overlap = max(0, min(cand_end, k_end) - max(cand_start, k_start))
        reciprocal = min(overlap / svlen, overlap / known_len)
        if reciprocal < MIN_RECIPROCAL:
            continue

        size_ratio = svlen / known_len
        if not (1 - SIZE_TOL) <= size_ratio <= (1 + SIZE_TOL):
            continue

        if reciprocal > best_score:
            best_score, best_name, best_itha = reciprocal, e["name"], e["ithaID"]

    if best_name:
        return "%s (%d%% reciprocal; ithaID=%s)" % (best_name, round(best_score * 100), best_itha)
    return "unknown_%s_%dbp" % (svtype, svlen)

