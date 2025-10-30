#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import argparse, math, sys, re
from collections import defaultdict, OrderedDict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from pycirclize import Circos

# =========================
# Palettes & simple helpers
# =========================

DEFAULT_KARYO_EDGE   = "#000000"
DEFAULT_KARYO_FILL   = "#808080"
DEFAULT_KARYO_PALETTE_NAME  = "tableau20"
DEFAULT_TRACK_PALETTE_NAME  = "selected"

PALETTES = {
    "tableau20":  ["#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7",
                   "#9C755F","#BAB0AC","#AEC7E8","#FFBB78","#FF9896","#98DF8A","#C5B0D5","#F7B6D2",
                   "#C49C94","#C7C7C7","#DBDB8D","#17BECF"],
    "selected":   ["#ed1b34","#946eb7","#7fd5ff","#00aaad","#36802d","#4e88c7"],
}

def get_palette(name: str, n: int) -> List[str]:
    base = PALETTES.get(name, [])
    if not base:
        return [DEFAULT_KARYO_FILL]*n
    return (base * (n // len(base) + 1))[:n]

def parse_list(s: Optional[str]) -> List[str]:
    return [t.strip() for t in s.split(",")] if s else []

def format_bp(x: int) -> str:
    if x >= 1_000_000_000: return f"{x/1e9:.0f} Gb"
    if x >= 1_000_000:     return f"{x/1e6:.0f} Mb"
    if x >= 1_000:         return f"{x/1e3:.0f} kb"
    return f"{x} bp"

def nice_ticks(seq_len: int, max_labels: int) -> Tuple[List[int], List[str]]:
    if seq_len <= 0: return [], []
    raw = max(1, seq_len/max_labels)
    exp = int(math.floor(math.log10(raw)))
    base = 10**exp
    step = next((s*base for s in (1,2,5,10) if s*base >= raw), 10*base)
    pos = list(range(0, seq_len+1, int(step))) or [0, seq_len]
    labels = [format_bp(p) for p in pos]
    if labels and pos[0] == 0:
        labels[0] = ""
    return pos, labels

# =========================
# Parsing / Binning helpers
# =========================

def fasta_lengths(fa: str) -> Dict[str,int]:
    lens, name, L = {}, None, 0
    with open(fa) as f:
        for ln in f:
            if ln.startswith(">"):
                if name is not None: lens[name] = L
                name, L = ln[1:].split()[0], 0
            else:
                L += len(ln.strip())
    if name is not None: lens[name] = L
    return lens

def gc_windows_from_fasta(fa: str, window: int, subset: Optional[set]) -> Dict[str, List[Tuple[int,int,float]]]:
    """
    Return {seq_name: [(start,end,gc%), ...]} with fixed-size windows.
    Combined version (no generator, no nested defs).
    """
    out: Dict[str, List[Tuple[int,int,float]]] = {}
    name: Optional[str] = None
    L = 0
    gc: List[int] = []
    tot: List[int] = []
    with open(fa) as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            if s.startswith(">"):
                # finalize previous sequence
                if name is not None and (subset is None or name in subset):
                    bins: List[Tuple[int,int,float]] = []
                    for i in range(len(tot)):
                        a = i*window + 1
                        b = min((i+1)*window, L)
                        pct = (100.0*gc[i]/tot[i]) if tot[i] else 0.0
                        bins.append((a,b,pct))
                    out[name] = bins
                # start new sequence
                name = s[1:].split()[0]
                L = 0; gc = []; tot = []
            else:
                for ch in s.upper():
                    L += 1
                    i = (L-1)//window
                    if i == len(tot):
                        gc.append(0); tot.append(0)
                    if ch in "ACGT":
                        tot[i]+=1
                        if ch in "GC": gc[i]+=1
    # finalize last sequence
    if name is not None and (subset is None or name in subset):
        bins = []
        for i in range(len(tot)):
            a = i*window + 1
            b = min((i+1)*window, L)
            pct = (100.0*gc[i]/tot[i]) if tot[i] else 0.0
            bins.append((a,b,pct))
        out[name] = bins
    return out

def read_gff_starts(
    gff: str,
    include_biotype: Optional[str],
    include_feature: Optional[str],
    exclude_biotype: Optional[str],
    subset: Optional[set]
) -> Dict[str, List[int]]:
    starts: Dict[str, List[int]] = defaultdict(list)
    with open(gff) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9: continue
            seq, feat, start_s = c[0], c[2], c[3]
            if subset and seq not in subset: continue
            if include_feature and feat != include_feature: continue
            try: start = int(start_s)
            except: continue
            attrs = dict(x.split("=",1) for x in c[8].split(";") if "=" in x)
            bt = attrs.get("gene_biotype") or attrs.get("biotype")
            if include_biotype and bt != include_biotype: continue
            if exclude_biotype and bt == exclude_biotype: continue
            starts[seq].append(start)
    return starts

def bins_from_starts(starts: List[int], window: int, seq_len: int) -> List[Tuple[int,int,int]]:
    n = math.ceil(seq_len / window)
    arr = [0]*n
    for s in starts:
        i = (s-1)//window
        if 0 <= i < n: arr[i]+=1
    return [(i*window+1, min((i+1)*window, seq_len), c) for i,c in enumerate(arr)]

def read_repeatmasker_multi(
    path: Optional[str], classes: Sequence[str], window: int, subset: Optional[set]
) -> Dict[str, Dict[str, List[Tuple[int,int,int]]]]:
    if not path: return {}
    acc: Dict[str, Dict[str, Dict[int, set]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("SW") or s.startswith("score"): continue
            flds = s.split()
            if len(flds) < 15: continue
            seq, repclass, rid = flds[4], flds[10], flds[14]
            if subset and seq not in subset: continue
            try: begin = int(flds[5])
            except: continue
            i = (begin-1)//window
            for cls in classes:
                if repclass.startswith(cls): acc[cls][seq][i].add(rid)
    out: Dict[str, Dict[str, List[Tuple[int,int,int]]]] = {}
    for cls, by_seq in acc.items():
        out[cls] = { seq: [(i*window+1,(i+1)*window,len(ids)) for i,ids in sorted(idx.items())]
                     for seq, idx in by_seq.items() }
    return out

# ============== RBH → link map ==============

def load_pairs_tsv(path: str) -> List[Tuple[str,str]]:
    pairs=[]
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            cols = ln.rstrip("\n").split("\t")
            if len(cols) >= 2: pairs.append((cols[0], cols[1]))
    return pairs

def load_map(path: Optional[str]) -> Dict[str,str]:
    if not path: return {}
    m={}
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            a,b = ln.rstrip("\n").split("\t",1)
            m[a]=b
    return m

def index_gff_midpoints(gff: str, subset: Optional[set]) -> Dict[str, Tuple[str,int,str]]:
    """Feature ID (mRNA/transcript/gene) and CDS Parent -> (seqid, midpoint, strand)."""
    idx={}
    with open(gff) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            c = ln.rstrip("\n").split("\t")
            if len(c)<9: continue
            seq, typ, s, e, strand = c[0], c[2], c[3], c[4], c[6]
            if subset and seq not in subset: continue
            try: s=int(s); e=int(e)
            except: continue
            attrs = dict(x.split("=",1) for x in c[8].split(";") if "=" in x)
            fid = attrs.get("ID")
            if typ in ("mRNA","transcript","gene") and fid:
                idx[fid]=(seq,(s+e)//2,strand)
            if typ=="CDS":
                pid=attrs.get("Parent")
                if pid and pid not in idx:
                    idx[pid]=(seq,(s+e)//2,strand)
    return idx

def auto_transform_candidates(pid: str) -> List[str]:
    _PT2MRNA_RE = re.compile(r"(.*)\.p([0-9]+)$")
    cands = [pid]
    m = _PT2MRNA_RE.match(pid)
    if m:
        cands.append(f"{m.group(1)}.t{m.group(2)}")
    if "." in pid:
        cands.append(pid.split(".",1)[0])
    out=[]; seen=set()
    for x in cands:
        if x not in seen:
            seen.add(x); out.append(x)
    return out

class Mapper:
    """Maps protein IDs to GFF feature IDs using explicit map and fallbacks (no closures)."""
    def __init__(self, explicit_map: Dict[str,str], gff_index: Dict[str, Tuple[str,int,str]]):
        self.explicit_map = explicit_map or {}
        self.gff_ids = set(gff_index.keys())

    def map_one(self, prot: str) -> Optional[str]:
        if prot in self.explicit_map:
            v = self.explicit_map[prot]
            return v if v in self.gff_ids else None
        for cand in auto_transform_candidates(prot):
            if cand in self.gff_ids:
                return cand
        return None

# =====================
# Labels / utilities
# =====================

def canonical_pair_key(seqid: str) -> str:
    """Normalize to 'chr...' (drop hap prefixes etc.) for color/label pairing."""
    s = seqid.lower()
    s = re.sub(r'^(hap(?:[ab]|[12])[_-])', '', s)
    m = re.search(r'(chr[_-]?\w+)$', s)
    return m.group(1) if m else s

def short_label(seqid: str) -> str:
    lab = canonical_pair_key(seqid)
    return lab.replace("chr_", "chr").replace("chr-", "chr")

def read_ids(path: Optional[str]) -> Optional[List[str]]:
    if not path: return None
    with open(path) as f:
        vals=[ln.strip() for ln in f if ln.strip()]
    return vals

# =============================================================
# Class: DualGenomeCircos (no nested defs)
# =============================================================

class DualGenomeCircos:
    def __init__(self, *, sector_space: float, group_gap_deg: float,
                 group_band_width: float, group_label_mode: str, group_label_pos: float,
                 group_label_offset: float, group_label_size: float,
                 outer_radius: float, karyo_thick: float, karyo_gap: float, track_gap: float,
                 xtick_interval: Optional[int], max_ticks_per_sector: int, tick_size: float, seq_label_size: float,
                 karyo_palette: str, track_palette: str,
                 groupA_color: str, groupB_color: str,
                 link_alpha: float):
        self.sector_space = sector_space
        self.group_gap_deg = group_gap_deg
        self.group_band_width = group_band_width
        self.group_label_mode = group_label_mode
        self.group_label_pos = group_label_pos
        self.group_label_offset = group_label_offset
        self.group_label_size = group_label_size
        self.outer_radius = outer_radius
        self.karyo_thick = karyo_thick
        self.karyo_gap = karyo_gap
        self.track_gap = track_gap
        self.xtick_interval = xtick_interval
        self.max_ticks_per_sector = max_ticks_per_sector
        self.tick_size = tick_size
        self.seq_label_size = seq_label_size
        self.karyo_palette = karyo_palette
        self.track_palette = track_palette
        self.groupA_color = groupA_color
        self.groupB_color = groupB_color
        self.link_alpha = link_alpha

    @staticmethod
    def ring_radii(r1: float, r2: float, track_heights: List[float], track_gap: float) -> List[Tuple[float,float]]:
        avail = r2 - r1 - track_gap*(len(track_heights)-1)
        widths = [avail * t/sum(track_heights) for t in track_heights]
        r=[]; cur=r1
        for w in widths:
            r.append((cur, cur+w))
            cur += w + track_gap
        return r

    @staticmethod
    def zero_deg_pad(track_obj) -> None:
        '''
        This function is applied to pad for space for line and strokes. 
        '''
        for attr in ("deg_pad", "angle_pad", "theta_pad"):
            if hasattr(track_obj, attr):
                try: setattr(track_obj, attr, 0.0)
                except Exception: pass

    def axis_no_pad(self, track_obj, **kw) -> None:
        self.zero_deg_pad(track_obj)
        track_obj.axis(**kw)

    @staticmethod
    def deg_lim_for(circos: Circos, names: List[str]) -> Optional[Tuple[float,float]]:
        ss=[s for s in circos.sectors if s.name in names]
        if not ss: return None
        d=[s.deg_lim for s in ss]
        return (min(a for a,b in d), max(b for a,b in d))

    @staticmethod
    def compute_gap_bp(total_bp: int, gap_deg_each: float, sector_space_deg: float, n_sectors_with_gaps: int
                       ) -> Tuple[int, float]:
        '''
        Purpose:
            Convert a requested angular gap size (in degrees) into a pseudo-sector
            length (in bp), so that the transparent gap wedges render with exactly
            the specified angular width in the circos plot. Also adjusts per-sector
            spacing so total layout fits into 360°.

        Variables:
            total_bp        - total bases from real chromosomes
            gap_deg_each    - angular size (deg) for one gap wedge (we place 2 gaps)
            sector_space_deg- desired spacing (deg) between every sector boundary
            n_sectors_with_gaps - total number of sectors including the 2 gaps

        Math:
            scale = (360 - total_space) / (total_bp + 2*G)
            require: G*scale = gap_deg_each
            solve for G → G = gap_deg_each * total_bp / (360 - total_space - 2*gap_deg_each)

        Returns:
            (gap_bp, adjusted_space_deg)
            gap_bp: pseudo length (bp) assigned to each gap sector
            adjusted_space_deg: corrected spacing per boundary if original spacing
                            would make the circle exceed 360°
        '''
        T = 360.0
        S = sector_space_deg * n_sectors_with_gaps
        denom = T - S - 2.0*gap_deg_each
        if denom <= 0:
            max_space = max(0.0, (T - 2.0*gap_deg_each) / n_sectors_with_gaps - 1e-6)
            S = max_space * n_sectors_with_gaps
            denom = T - S - 2.0*gap_deg_each
        G = (gap_deg_each * total_bp) / max(1e-6, denom)
        return max(1, int(round(G))), S / n_sectors_with_gaps

    def build_sectors_with_gaps(self, A_sizes: Dict[str,int], B_sizes: Dict[str,int]
                                ) -> Tuple[OrderedDict, Dict[str,bool], str, str, float]:
        A_names = list(A_sizes.keys()); B_names = list(B_sizes.keys())
        total_bp_real = sum(A_sizes.values()) + sum(B_sizes.values())
        N = len(A_names) + len(B_names) + 2  # +2 gap sectors
        gap_bp, sector_space = self.compute_gap_bp(total_bp_real, self.group_gap_deg, self.sector_space, N)
        GAP1 = "__GAP__A_to_B__"; GAP2 = "__GAP__B_to_A__"
        sectors = OrderedDict()
        for k in A_names: sectors[k] = A_sizes[k]
        sectors[GAP1] = gap_bp
        for k in reversed(B_names): sectors[k] = B_sizes[k]
        sectors[GAP2] = gap_bp
        sector2clockwise = {k: True for k in A_names}
        sector2clockwise.update({k: False for k in B_names})
        sector2clockwise[GAP1] = True; sector2clockwise[GAP2] = True
        return sectors, sector2clockwise, GAP1, GAP2, sector_space

    def draw_karyotype_and_label(self, sector, k1, k2, color, label):
        '''
        k1 = inner radius; k2 = outer radius
        label text on the karyotype band.
        rotation=mid+90 to make text orthogonal to radius.
        '''
        kt = sector.add_track((k1, k2))
        self.axis_no_pad(kt, fc=color, ec=DEFAULT_KARYO_EDGE)
        mid = 0.5*(sector.deg_lim[0]+sector.deg_lim[1])
        kt.text(label, r=(k1+k2)/2, size=self.seq_label_size, fontweight="bold", fontstyle="italic",
                rotation=mid+90, rotation_mode="anchor", va="center", ha="center", zorder=50)

    def draw_tick_ring(self, sector, tick_r1, tick_r2, L: int):
        '''
        zero_deg_pad to avoid tick cropping/addition at 0°/360° boundary.
        L = the sequence length (in base pairs) for this sector.
        self.max_ticks_per_sector = the maximum number of labeled ticks you want.
        '''
        tkt = sector.add_track((tick_r1, tick_r2))
        self.zero_deg_pad(tkt)
        if self.xtick_interval:
            tkt.xticks_by_interval(self.xtick_interval, label_size=self.tick_size,
                                   label_orientation="vertical",
                                   label_formatter=format_bp, tick_length=0.9)
        else:
            pos, lab = nice_ticks(L, self.max_ticks_per_sector)
            if pos:
                tkt.xticks(pos, labels=lab, label_size=self.tick_size,
                           label_orientation="vertical", tick_length=0.9)

    # ============================
    # MODIFIED: global, uniform y-scales for genes & repeats (no percentiles)
    # ============================
    def draw_data_tracks(self, sector, rings, tracks, track_colors, gc_bins, gene_bins, rep_bins, seq_len):
        '''
        Draw per-sector data tracks (GC%, genes, RepeatMasker classes) using
        CONSISTENT y-scales across BOTH haplotypes and ALL chromosomes.

        Scaling rules:
          - GC: fixed 0–100.
          - Genes: global cap across BOTH haplotypes = max(count) over ALL windows of ALL chromosomes.
          - Repeat classes: per-class global cap across BOTH haplotypes = max(count) over ALL windows of ALL chromosomes.
          - No percentile scaling.
        '''
        # Pull precomputed global caps (set in plot()) with safe fallbacks.
        global_gene_ymax = getattr(self, "global_gene_ymax", 1)
        if not isinstance(global_gene_ymax, (int, float)) or global_gene_ymax <= 0:
            global_gene_ymax = 1.0

        global_rep_caps = getattr(self, "global_rep_caps", {})
        if not isinstance(global_rep_caps, dict):
            global_rep_caps = {}

        for (token, (r1, r2)), color in zip(zip(tracks, rings), track_colors):
            t = sector.add_track((r1, r2), r_pad_ratio=0.06)
            self.axis_no_pad(t)
            t.axis()

            if token == "gc":
                bins = gc_bins or []
                ymax = 100.0

            elif token == "genes":
                bins = gene_bins or []
                ymax = float(global_gene_ymax)

            elif token.startswith("repeat:"):
                cls = token.split(":", 1)[1]
                binmap = rep_bins.get(cls, {}) if rep_bins else {}   # {seqid: [(a,b,c), ...]}
                bins   = binmap.get(sector.name, [])                 # this sector's bins
                cap = global_rep_caps.get(cls, 1)
                if not isinstance(cap, (int, float)) or cap <= 0: #isinstance(3, (int, float)) True
                    cap = 1.0
                ymax = float(cap)

            else:
                bins = []
                ymax = 1.0

            if bins:
                # Clip to seq_len; compute centers and inclusive widths
                xs = np.array([0.5 * (a + min(b, seq_len)) for a, b, _ in bins if a <= seq_len])
                ys = np.array([v for a, b, v in bins if a <= seq_len])
                ws = np.array([min(b, seq_len) - a + 1 for a, b, _ in bins if a <= seq_len])

                if xs.size:
                    t.bar(xs, np.clip(ys, 0, ymax), width=ws, ec="none", fc=color)

    def draw_group_segment(self, circos: Circos, k2: float, label: str,
                           deg_lim: Tuple[float,float], color: str) -> None:
        band_inner = k2
        band_outer = k2 + self.group_band_width
        circos.rect(r_lim=(band_inner, band_outer), deg_lim=deg_lim,
                    fc=color, ec="black", lw=0.6, alpha=0.9, zorder=20)
        center = 0.5*(deg_lim[0]+deg_lim[1])
        if self.group_label_mode == "onband":
            r = band_inner + self.group_label_pos*(band_outer-band_inner) + self.group_label_offset
            circos.text(label, r=r, deg=center, adjust_rotation=True,
                        size=self.seq_label_size, fontweight="bold", fontstyle="italic",
                        va="center", ha="center", zorder=120, clip_on=False)
        else:
            tick_r1 = k2 + self.group_band_width
            tick_r2 = tick_r1 + 0.9
            label_r = tick_r2 + 3.5 + self.group_label_offset
            circos.text(label, r=label_r, deg=center, adjust_rotation=True,
                        size=self.group_label_size, fontweight="bold",
                        va="center", ha="center", zorder=120, clip_on=False)

    def plot(self, *, out_png: str,
             A_sizes: Dict[str,int], A_gc, A_gene_bins, A_rep_bins,
             B_sizes: Dict[str,int], B_gc, B_gene_bins, B_rep_bins,
             tracks: List[str], track_heights: List[float],
             links: Optional[List[Tuple[str,int,str,int,bool]]]) -> None:
        # ============================
        # MODIFIED: precompute GLOBAL caps before any drawing
        # ============================
        # Global gene ymax across BOTH haplotypes
        '''
        *: from now on, all parameters must be passed by name.

        '''
        global_gene_max = 0
        for bins_by_seq in (A_gene_bins, B_gene_bins):
            for arr in bins_by_seq.values():           # arr: [(a,b,count), ...]
                for _a, _b, c in arr:
                    if c > global_gene_max:
                        global_gene_max = c
        self.global_gene_ymax = global_gene_max if global_gene_max > 0 else 1

        # Global repeat caps per class across BOTH haplotypes
        rep_caps: Dict[str, int] = {}
        for rep_dict in (A_rep_bins, B_rep_bins):      # each: {class: {seqid: [(a,b,count), ...]}}
            for cls, byseq in rep_dict.items():
                for arr in byseq.values():
                    for _a, _b, c in arr:
                        if c > rep_caps.get(cls, 0):
                            rep_caps[cls] = c
        # ensure >=1 to avoid degenerate scales
        self.global_rep_caps = {cls: (v if v > 0 else 1) for cls, v in rep_caps.items()}

        sectors, sector2clockwise, GAP1, GAP2, sector_space = self.build_sectors_with_gaps(A_sizes, B_sizes)
        circos = Circos(sectors, start=0.0, end=360.0, space=sector_space, endspace=True,
                        sector2clockwise=sector2clockwise)

        tick_ring_width = 0.9
        k2 = self.outer_radius - (self.group_band_width + tick_ring_width)
        k1 = k2 - self.karyo_thick
        data_outer = k1 - self.karyo_gap
        total_th = sum(track_heights) + self.track_gap*(len(track_heights)-1)
        data_r1, data_r2 = data_outer - total_th, data_outer
        rings = self.ring_radii(data_r1, data_r2, track_heights, self.track_gap)

        A_names = list(A_sizes.keys()); B_names = list(B_sizes.keys())
        pair_keys=[]
        for nm in A_names + B_names:
            pk = canonical_pair_key(nm)
            if pk not in pair_keys: pair_keys.append(pk)
        pair_colors = get_palette(self.karyo_palette, len(pair_keys))
        key2color = {k: pair_colors[i] for i,k in enumerate(pair_keys)}
        track_colors = get_palette(self.track_palette, len(tracks))

        gap_names = {GAP1, GAP2}
        for sector in circos.sectors:
            seq, L = sector.name, sector.size
            if seq in gap_names:
                continue

            self.draw_karyotype_and_label(sector, k1, k2, key2color[canonical_pair_key(seq)], short_label(seq))

            tick_r1 = k2 + self.group_band_width
            tick_r2 = tick_r1 + tick_ring_width
            self.draw_tick_ring(sector, tick_r1, tick_r2, L)

            if seq in A_sizes:
                gc_bins = A_gc.get(seq, []); gene_bins = A_gene_bins.get(seq, []); rep_bins = A_rep_bins; seq_len=A_sizes[seq]
            else:
                gc_bins = B_gc.get(seq, []); gene_bins = B_gene_bins.get(seq, []); rep_bins = B_rep_bins; seq_len=B_sizes[seq]

            self.draw_data_tracks(sector, rings, tracks, track_colors, gc_bins, gene_bins, rep_bins, seq_len)

        hap1_lim = self.deg_lim_for(circos, A_names)
        hap2_lim = self.deg_lim_for(circos, B_names)
        if hap1_lim: self.draw_group_segment(circos, k2, "hap1", hap1_lim, self.groupA_color)
        if hap2_lim: self.draw_group_segment(circos, k2, "hap2", hap2_lim, self.groupB_color)

        if links:
            link_r = max(0.0, data_r1 - 2.0)
            for s1,p1,s2,p2,_inv in links:
                color = key2color.get(canonical_pair_key(s1), "#888888")
                circos.link((s1,p1,p1), (s2,p2,p2), color=color, r1=link_r, r2=link_r, alpha=self.link_alpha)

        fig = circos.plotfig()
        fig.savefig(out_png, dpi=450, bbox_inches="tight")
        print(f"[OK] wrote {out_png}")

# =========================
# CLI wrapper
# =========================

def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Two-genome circos with transparent pseudo-sectors as group gaps and on-band group labels (class-based, no nested defs)."
    )
    # Genome A (outer: hap1)
    ap.add_argument("--A_fasta", required=True)
    ap.add_argument("--A_gff",   required=True)
    ap.add_argument("--A_repeatout", help="RepeatMasker .out for A")
    # Genome B (inner: hap2)
    ap.add_argument("--B_fasta", required=True)
    ap.add_argument("--B_gff",   required=True)
    ap.add_argument("--B_repeatout", help="RepeatMasker .out for B")

    # Tracks & data
    ap.add_argument("--tracks", required=True, help="Comma list INNER→OUTER, e.g. 'gc,repeat:LINE,genes'")
    ap.add_argument("--track_thicknesses", required=True, help="Comma floats (same count as --tracks)")
    ap.add_argument("--window", type=int, default=100000)
    ap.add_argument("--feature", help="GFF feature type to count (e.g., gene or mRNA)")
    ap.add_argument("--biotype", help="Only include features with gene_biotype==this value")
    ap.add_argument("--exclude_biotype", help="Exclude features with gene_biotype==this value")

    # SeqID filters
    ap.add_argument("--seqids", help="Apply this ID list to BOTH haplotypes unless overridden below")
    ap.add_argument("--A_seqids", help="Optional list (one per line) to subset/reorder A")
    ap.add_argument("--B_seqids", help="Optional list (one per line) to subset/reorder B")

    # RBH links
    ap.add_argument("--rbh_tsv", required=True, help="MMseqs2 RBH table (first two columns query,target)")
    ap.add_argument("--A_map", help="2-col TSV mapping: hap1_protein_id -> GFF feature ID")
    ap.add_argument("--B_map", help="2-col TSV mapping: hap2_protein_id -> GFF feature ID")
    ap.add_argument("--max_links", type=int, default=50000)
    ap.add_argument("--link_alpha", type=float, default=0.4)

    # Plot styling (full circle; no start/end angles)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--sector_space", type=float, default=0.2,
                    help="Small constant spacing (degrees) between EVERY adjacent sector")
    ap.add_argument("--karyotype_thickness", type=float, default=4.0)
    ap.add_argument("--karyotype_gap", type=float, default=1.5)
    ap.add_argument("--track_gap", type=float, default=2.0)
    ap.add_argument("--xtick_interval", type=int)
    ap.add_argument("--max_ticks_per_sector", type=int, default=6)
    ap.add_argument("--tick_label_size", type=float, default=5.5)
    ap.add_argument("--karyo_palette", default=DEFAULT_KARYO_PALETTE_NAME)
    ap.add_argument("--track_palette", default=DEFAULT_TRACK_PALETTE_NAME)
    ap.add_argument("--seq_label_size", type=float, default=6.0)

    # Group appearance
    ap.add_argument("--groupA_color", default="skyblue")
    ap.add_argument("--groupB_color", default="salmon")
    ap.add_argument("--group_band_width", type=float, default=2.5)
    ap.add_argument("--group_label_mode", choices=["onband","outside"], default="onband")
    ap.add_argument("--group_label_pos", type=float, default=0.55,
                    help="0=band inner edge, 1=band outer edge (lower = more inward)")
    ap.add_argument("--group_label_offset", type=float, default=0.0)
    ap.add_argument("--group_label_size", type=float, default=13.0)

    # Group gap (REAL pseudo-sectors)
    ap.add_argument("--group_gap_deg", type=float, default=12.0,
                    help="Degrees for EACH transparent gap sector between hap1↔hap2 and at the wrap")

    args = ap.parse_args()

    both_ids = read_ids(args.seqids)
    A_order = read_ids(args.A_seqids) or both_ids
    B_order = read_ids(args.B_seqids) or both_ids

    tracks = [t for t in parse_list(args.tracks) if t]
    thicks = [float(x) for x in parse_list(args.track_thicknesses)]
    if len(thicks) != len(tracks):
        sys.exit("Error: --track_thicknesses count must match --tracks")

    print("[INFO] Reading FASTA lengths…")
    A_all = fasta_lengths(args.A_fasta)
    B_all = fasta_lengths(args.B_fasta)
    A_names = A_order if A_order else sorted(A_all)
    B_names = B_order if B_order else sorted(B_all)
    A_sizes = OrderedDict((x, A_all[x]) for x in A_names if x in A_all)
    B_sizes = OrderedDict((x, B_all[x]) for x in B_names if x in B_all)
    A_set, B_set = set(A_sizes), set(B_sizes)

    print("[INFO] Computing GC% windows…")
    W = args.window
    A_gc = gc_windows_from_fasta(args.A_fasta, W, A_set)
    B_gc = gc_windows_from_fasta(args.B_fasta, W, B_set)

    print("[INFO] Binning gene starts…")
    A_starts = read_gff_starts(args.A_gff, args.biotype, args.feature, args.exclude_biotype, A_set)
    B_starts = read_gff_starts(args.B_gff, args.biotype, args.feature, args.exclude_biotype, B_set)
    A_gene_bins = {seq: bins_from_starts(A_starts.get(seq,[]), W, A_sizes[seq]) for seq in A_sizes}
    B_gene_bins = {seq: bins_from_starts(B_starts.get(seq,[]), W, B_sizes[seq]) for seq in B_sizes}

    print("[INFO] Reading RepeatMasker…")
    rep_classes = [t.split(":",1)[1] for t in tracks if t.startswith("repeat:")]
    A_rep_bins = read_repeatmasker_multi(args.A_repeatout, rep_classes, W, A_set)
    B_rep_bins = read_repeatmasker_multi(args.B_repeatout, rep_classes, W, B_set)

    print("[INFO] Building RBH link list…")
    pairs = load_pairs_tsv(args.rbh_tsv)[:args.max_links]
    A_idx = index_gff_midpoints(args.A_gff, A_set)
    B_idx = index_gff_midpoints(args.B_gff, B_set)
    A_map = load_map(args.A_map)
    B_map = load_map(args.B_map)

    A_mapper = Mapper(A_map, A_idx)
    B_mapper = Mapper(B_map, B_idx)

    links=[]
    skipped=0
    for q,t in pairs:
        fq = A_mapper.map_one(q); ft = B_mapper.map_one(t)
        if fq is None or ft is None or fq not in A_idx or ft not in B_idx:
            skipped+=1; continue
        s1,p1,str1 = A_idx[fq]
        s2,p2,str2 = B_idx[ft]
        if s1 not in A_set or s2 not in B_set: continue
        inv = (str1 != str2)
        links.append((s1,p1,s2,p2,inv))
    if skipped:
        print(f"[WARN] RBH pairs without GFF mapping skipped: {skipped}", file=sys.stderr)
    print(f"[INFO] RBH links retained: {len(links)}")

    plotter = DualGenomeCircos(
        sector_space=args.sector_space, group_gap_deg=args.group_gap_deg,
        group_band_width=args.group_band_width, group_label_mode=args.group_label_mode,
        group_label_pos=args.group_label_pos, group_label_offset=args.group_label_offset,
        group_label_size=args.group_label_size, outer_radius=100.0,
        karyo_thick=args.karyotype_thickness, karyo_gap=args.karyotype_gap, track_gap=args.track_gap,
        xtick_interval=args.xtick_interval, max_ticks_per_sector=args.max_ticks_per_sector,
        tick_size=args.tick_label_size, seq_label_size=args.seq_label_size,
        karyo_palette=args.karyo_palette, track_palette=args.track_palette,
        groupA_color=args.groupA_color, groupB_color=args.groupB_color,
        link_alpha=args.link_alpha,
    )

    plotter.plot(
        out_png=args.out_png,
        A_sizes=A_sizes, A_gc=A_gc, A_gene_bins=A_gene_bins, A_rep_bins=A_rep_bins,
        B_sizes=B_sizes, B_gc=B_gc, B_gene_bins=B_gene_bins, B_rep_bins=B_rep_bins,
        tracks=tracks, track_heights=thicks, links=links,
    )

if __name__ == "__main__":
    main()
