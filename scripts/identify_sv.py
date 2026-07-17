KNOWN_SVS = {
    # name: chrom, expected_start, expected_end, type
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
    pos = int(pos)
    svlen = abs(int(svlen)) if svlen else 0
    sv_end = pos + svlen
    
    best_match = None
    best_overlap = 0
    
    for name, info in KNOWN_SVS.items():
        if chrom != info["chrom"]:
            continue
        if svtype != info["type"]:
            continue
        
        # Calculate overlap between detected SV and known SV
        overlap_start = max(pos, info["del_start"])
        overlap_end = min(sv_end, info["del_end"])
        overlap = max(0, overlap_end - overlap_start)
        
        # What fraction of the detected SV overlaps the known SV?
        if svlen > 0:
            overlap_fraction = overlap / svlen
        else:
            overlap_fraction = 0
        
        if overlap_fraction > 0.5 and overlap_fraction > best_overlap:
            best_overlap = overlap_fraction
            best_match = name
    
    if best_match:
        return f"{best_match} ({round(best_overlap * 100)}% overlap)"
    return f"unknown_{svtype}_{svlen}bp"
