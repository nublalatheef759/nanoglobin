#!/usr/bin/env python3
"""
copy-number interpretation for the HBA/HBB coverage ratio.

Maps the normalised depth ratio to an estimated copy number. Losses (deletions)
are crisp because the steps are large (1.0 -> 0.5 -> 0.0). Gains are subtler:
a het triplication is 3 copies vs a normal 2, i.e. ratio ~1.5 vs 1.0 -- a smaller
step than a deletion -- so gain bands need wider tolerance and are reported with
lower confidence than losses.

Normal diploid = 2 copies = ratio 1.0. Each copy is 0.5 of the diploid depth.
"""

def call_copy_number(ratio):
    """Return (copies, call, confidence) for a normalised depth ratio."""
    if ratio is None or ratio != ratio:   # NaN
        return (None, "undetermined", "none")

    # bands centred on n*0.5 (n copies). Loss bands are tight; gain bands wider.
    bands = [
        (0.00, 0.25, 0, "homozygous deletion (0 copies)",   "high"),
        (0.25, 0.75, 1, "heterozygous deletion (1 copy)",   "high"),
        (0.75, 1.25, 2, "normal (2 copies)",                "high"),
        (1.25, 1.75, 3, "triplication (3 copies)",          "moderate"),
        (1.75, 2.50, 4, "quadruplication (4 copies)",       "moderate"),
    ]
    for lo, hi, n, call, conf in bands:
        if lo <= ratio < hi:
            return (n, call, conf)
    if ratio >= 2.50:
        return (5, "high-level amplification (>=5 copies)", "low")
    return (None, "undetermined", "none")


if __name__ == "__main__":
    # test across the full range, both losses and gains
    tests = [
        (0.02, 0, "homozygous deletion"),
        (0.50, 1, "heterozygous deletion"),
        (0.87, 2, "normal"),           # real HET_a37 whole-window mean was 0.87
        (1.01, 2, "normal"),
        (1.50, 3, "triplication"),     # het aaa/aa
        (1.48, 3, "triplication"),
        (2.00, 4, "quadruplication"),
        (0.00, 0, "hom del"),
    ]
    print("%-6s %-6s %-32s %s" % ("ratio", "copies", "call", "confidence"))
    ok = True
    for ratio, exp_n, _ in tests:
        n, call, conf = call_copy_number(ratio)
        flag = "" if n == exp_n else "  <-- MISMATCH"
        if n != exp_n: ok = False
        print("%-6.2f %-6s %-32s %s%s" % (ratio, n, call, conf, flag))
    print("\nALL PASS" if ok else "SOME FAILED")
    