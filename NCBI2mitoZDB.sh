# usage: bash convert_one_line.sh INPUT.faa "Genus_species" > OUTPUT_for_MitoZ.faa
python3 - <<'PY' "$1" "$2"
import sys, re
inp = sys.argv[1]
species = (sys.argv[2] if len(sys.argv) > 2 else "Genus_sp.").replace(" ", "_")

gene_map = {
    "cox1":"COX1","cox2":"COX2","cox3":"COX3","cob":"CYTB",
    "nad1":"ND1","nad2":"ND2","nad3":"ND3","nad4":"ND4","nad4l":"ND4L","nad5":"ND5","nad6":"ND6",
    "atp6":"ATP6","atp8":"ATP8"
}

def flush(h, seq):
    if not h: return
    s = "".join(seq).replace("\n","").replace(" ","")
    L = len(s)
    mg = re.search(r"\[gene=([^\]]+)\]", h, re.I)
    g = mg.group(1).lower() if mg else "UNKNOWN"
    gene = gene_map.get(g, g.upper())
    mp = re.search(r"\[protein_id=([^\]]+)\]", h)
    if mp:
        acc = re.sub(r"\..*$", "", mp.group(1))
    else:
        m2 = re.search(r"_prot_([A-Za-z0-9_.-]+)", h)
        acc = re.sub(r"\..*$", "", m2.group(1)) if m2 else "UNKNOWN"
    print(f">gi_NC_{acc}_{gene}_{species}_{L}_aa")
    print(s)  # one line only

header, seq = None, []
with open(inp) as f:
    for line in f:
        if line.startswith(">"):
            flush(header, seq)
            header, seq = line.strip(), []
        else:
            seq.append(line.strip())
flush(header, seq)
PY
