"""
SV identification against a catalogue of known thalassaemia structural variants.

A candidate SV is only assigned a name if it satisfies BOTH:
  1. Reciprocal overlap >= MIN_RECIPROCAL in both directions, i.e. the candidate
     lies mostly within the known variant AND the known variant lies mostly within
     the candidate. One-directional containment is not sufficient: a 92bp deletion
     sits entirely inside the Hb Lepore window but is not Hb Lepore.
  2. Size concordance: |candidate| / |known| within (1 +/- SIZE_TOL). This is what
     separates -a3.7 (~3.8kb) from -a4.2 (~4.3kb), whose windows overlap by 89%
     and which reciprocal overlap alone cannot distinguish.

Breakpoint-proximity filtering is deliberately NOT applied: the coordinates in
KNOWN_SVS are unsourced and disagree with published breakpoints by ~900bp (see
README). Adding a breakpoint filter on top of uncertain coordinates would convert
false positives into false negatives.
"""

MIN_RECIPROCAL = 0.5
SIZE_TOL = 0.10

KNOWN_SVS = {
    "-a3.7":      {"chrom": "chr16", "del_start": 172800, "del_end": 176600, "type": "DEL"},
    "-a4.2":      {"chrom": "chr16", "del_start": 173200, "del_end": 177500, "type": "DEL"},
    "--SEA":      {"chrom": "chr16", "del_start": 158000, "del_end": 178000, "type": "DEL"},
    "--FIL":      {"chrom": "chr16", "del_start": 144000, "del_end": 178000, "type": "DEL"},
    "--MED-I":    {"chrom": "chr16", "del_start": 155000, "del_end": 172500, "type": "DEL"},
    "--THAI":     {"chrom": "chr16", "del_start": 151000, "del_end": 178000, "type": "DEL"},
    "anti-3.7":   {"chrom": "chr16", "del_start": 172800, "del_end": 176600, "type": "DUP"},
    "anti-4.2":   {"chrom": "chr16", "del_start": 173200, "del_end": 177500, "type": "DUP"},
    "619bp del":  {"chrom": "chr11", "del_start": 5225464, "del_end": 5226100, "type": "DEL"},
    "Hb Lepore":  {"chrom": "chr11", "del_start": 5225400, "del_end": 5232000, "type": "DEL"},
}


def identify_sv(chrom, pos, svtype, svlen):
    """Return a known variant name, or unknown_<TYPE>_<LEN>bp if nothing qualifies."""
    pos = int(pos)
    svlen = abs(int(svlen)) if svlen else 0
    if svlen == 0:
        return "unknown_%s_0bp" % svtype
    sv_end = pos + svlen

    best_match = None
    best_score = 0.0

    for name, info in KNOWN_SVS.items():
        if chrom != info["chrom"] or svtype != info["type"]:
            continue

        known_len = info["del_end"] - info["del_start"]
        if known_len <= 0:
            continue

        overlap = max(0, min(sv_end, info["del_end"]) - max(pos, info["del_start"]))

        # reciprocal overlap: the weaker of the two directions must clear the bar
        reciprocal = min(overlap / svlen, overlap / known_len)
        if reciprocal < MIN_RECIPROCAL:
            continue

        # size concordance
        size_ratio = svlen / known_len
        if not (1 - SIZE_TOL) <= size_ratio <= (1 + SIZE_TOL):
            continue

        if reciprocal > best_score:
            best_score = reciprocal
            best_match = name

    if best_match:
        return "%s (%d%% reciprocal)" % (best_match, round(best_score * 100))
    return "unknown_%s_%dbp" % (svtype, svlen)

