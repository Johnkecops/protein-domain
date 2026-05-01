#!/usr/bin/env python3
"""
Protein Domain Analysis Toolkit
================================
Author  : Dr. Arli Aditya Parikesit
Affil.  : Department of Bioinformatics, i3L University, Jakarta
ORCID   : 0000-0001-8716-3926
Date    : 2026

Integrates all legacy shell scripts into a single, menu-driven Python
pipeline with a real S. cerevisiae use-case mined live from UniProt.

Tools integrated
----------------
1. Domain co-occurrence analysis       (Domain cooccurence.sh / -o.sh)
2. BED coordinate → FASTA extraction  (Convert BED coordinate to FASTA.sh)
3. Redundant FASTA removal             (Remove Redundant FASTA.pl)
4. Duplicate detection in bio-networks (Duplicates in Biological Network.sh)
5. Overlapping domain region collapse  (Domain cooccurence.sh core logic)
6. UniProt proteomics data mining      (new — REST API)
7. XMGrace .agr co-occurrence output   (new — publication-ready)

References
----------
Parikesit et al. (2014). Malaysian J. Fund. Appl. Sci., 10(2), 65-75.
  https://doi.org/10.11113/mjfas.v10n2.57
Parikesit et al. (2011). Genes, 2(4), 912-924.
  https://doi.org/10.3390/genes2040912
"""

from __future__ import annotations

import itertools
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns

# ── optional heavy dep ──────────────────────────────────────────────────────
try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

try:
    from pyfaidx import Fasta
    HAS_PYFAIDX = True
except ImportError:
    HAS_PYFAIDX = False

# ── constants ───────────────────────────────────────────────────────────────
UNIPROT_API      = "https://rest.uniprot.org/uniprotkb/search"
SCEREVISIAE_TAXID = 559292
RESULTS_DIR       = Path("results")


# ════════════════════════════════════════════════════════════════════════════
#  1. UniProt mining
# ════════════════════════════════════════════════════════════════════════════

def mine_uniprot(
    organism_id: int  = SCEREVISIAE_TAXID,
    reviewed:    bool = True,
    max_results: int  = 500,
    output_dir:  Path = RESULTS_DIR,
) -> pd.DataFrame:
    """
    Mine protein domain annotations from UniProt via REST API.

    Fields fetched: accession, protein name, gene names, organism,
    sequence length, domain features (ft_domain), GO biological process.

    Parameters
    ----------
    organism_id : NCBI Taxonomy ID  (559292 = S. cerevisiae)
    reviewed    : Swiss-Prot only when True
    max_results : upper bound on proteins returned
    output_dir  : where to write the TSV cache

    Returns
    -------
    pd.DataFrame  with one row per protein
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    query = f"organism_id:{organism_id}"
    if reviewed:
        query += " AND reviewed:true"

    fields = (
        "accession,protein_name,gene_names,organism_name,"
        "length,ft_domain,go_p"
    )

    params: dict | None = {
        "query":  query,
        "format": "tsv",
        "fields": fields,
        "size":   min(max_results, 500),
    }

    print(f"  [UniProt] Querying organism={organism_id}, reviewed={reviewed}, "
          f"max={max_results} ...")

    rows:   list[list[str]] = []
    header: list[str]       = []
    url = UNIPROT_API

    while url and len(rows) < max_results:
        resp = requests.get(url, params=params, timeout=60)
        resp.raise_for_status()

        lines = resp.text.strip().split("\n")
        if len(lines) < 2:
            break

        if not header:
            header = lines[0].split("\t")

        for line in lines[1:]:
            if len(rows) >= max_results:
                break
            rows.append(line.split("\t"))

        # pagination via Link header
        link = resp.headers.get("Link", "")
        next_url: str | None = None
        for part in link.split(","):
            if 'rel="next"' in part:
                m = re.search(r"<(.+?)>", part)
                if m:
                    next_url = m.group(1)
        url    = next_url
        params = None          # params only used on first request
        if next_url:
            time.sleep(0.3)

    df = pd.DataFrame(rows, columns=header)

    # Pad short rows (API occasionally returns ragged lines)
    for col in header:
        if col not in df.columns:
            df[col] = ""

    out = output_dir / f"uniprot_{organism_id}.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"  [UniProt] {len(df)} proteins saved → {out}")
    return df


# ════════════════════════════════════════════════════════════════════════════
#  2. Domain parsing + co-occurrence
# ════════════════════════════════════════════════════════════════════════════

def _parse_domains(raw: str) -> list[str]:
    """
    Extract domain names from a UniProt ft_domain cell.

    UniProt TSV domain format (structured):
        DOMAIN 1..100; /note="Kinase catalytic"; DOMAIN 110..200; /note="SH2"
    Fallback (semicolon list):
        Kinase catalytic; SH2; WD repeat
    """
    if not raw or (isinstance(raw, float) and np.isnan(raw)):
        return []
    raw = str(raw).strip()
    if not raw:
        return []

    notes = re.findall(r'/note="([^"]+)"', raw)
    if notes:
        return [n.strip() for n in notes if n.strip()]

    # fallback — plain semicolon list
    parts = [p.strip() for p in raw.split(";")]
    return [p for p in parts if p and not re.fullmatch(r"[\d\.\s]+", p)]


def compute_domain_cooccurrence(df: pd.DataFrame) -> tuple[pd.DataFrame, Counter]:
    """
    Build a symmetric domain co-occurrence matrix.

    For every protein that carries ≥2 annotated domains, all unordered
    pairs (d_i, d_j) are counted.  Self-pairs (d, d) record per-domain
    frequency.

    Returns
    -------
    matrix      : pd.DataFrame  (domains × domains, integer counts)
    cooccurrence: Counter keyed by (domain_a, domain_b) where a ≤ b
    """
    domain_col = next(
        (c for c in df.columns if "domain" in c.lower()),
        None,
    )
    if domain_col is None:
        raise ValueError(
            "No domain column found.  Expected a column whose name "
            "contains 'domain' (e.g. 'ft_domain')."
        )

    cooccurrence: Counter = Counter()

    for _, row in df.iterrows():
        domains = _parse_domains(row[domain_col])
        if not domains:
            continue
        unique = sorted(set(domains))

        # single-domain count (diagonal)
        for d in unique:
            cooccurrence[(d, d)] += 1

        # pairwise co-occurrence
        for d1, d2 in itertools.combinations(unique, 2):
            pair = (min(d1, d2), max(d1, d2))
            cooccurrence[pair] += 1

    all_domains = sorted({d for pair in cooccurrence for d in pair})
    if not all_domains:
        return pd.DataFrame(), cooccurrence

    matrix = pd.DataFrame(0, index=all_domains, columns=all_domains, dtype=int)
    for (d1, d2), cnt in cooccurrence.items():
        matrix.loc[d1, d2] = cnt
        matrix.loc[d2, d1] = cnt

    return matrix, cooccurrence


# ════════════════════════════════════════════════════════════════════════════
#  3. XMGrace .agr output
# ════════════════════════════════════════════════════════════════════════════

_AGR_FONT_MAP = """\
@map font 0 to "Times-Roman", "Times-Roman"
@map font 1 to "Times-Italic", "Times-Italic"
@map font 2 to "Times-Bold", "Times-Bold"
@map font 3 to "Times-BoldItalic", "Times-BoldItalic"
@map font 4 to "Helvetica", "Helvetica"
@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"
@map font 6 to "Helvetica-Bold", "Helvetica-Bold"
@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"
@map font 8 to "Symbol", "Symbol"
@map font 9 to "ZapfDingbats", "ZapfDingbats"
@map font 10 to "Courier", "Courier"
@map font 11 to "Courier-Oblique", "Courier-Oblique"
@map font 12 to "Courier-Bold", "Courier-Bold"
@map font 13 to "Courier-BoldOblique", "Courier-BoldOblique"
"""

_AGR_COLOR_MAP = """\
@map color 0 to (255, 255, 255), "white"
@map color 1 to (0, 0, 0), "black"
@map color 2 to (255, 0, 0), "red"
@map color 3 to (0, 200, 0), "green"
@map color 4 to (0, 0, 255), "blue"
@map color 5 to (255, 255, 0), "yellow"
@map color 6 to (188, 143, 143), "brown"
@map color 7 to (220, 220, 220), "grey"
@map color 8 to (148, 0, 211), "violet"
@map color 9 to (0, 255, 255), "cyan"
@map color 10 to (255, 127, 0), "orange"
@map color 11 to (64, 224, 208), "turquoise"
@map color 12 to (255, 0, 255), "magenta"
@map color 13 to (100, 149, 237), "indigo"
@map color 14 to (255, 105, 180), "maroon"
@map color 15 to (143, 143, 188), "periwinkle"
"""


def generate_xmgrace_cooccurrence(
    cooccurrence: Counter,
    output_path:  Path | str,
    top_n:        int = 20,
    title:        str = "Domain Co-occurrence",
) -> str:
    """
    Write an XMGrace project file (.agr) showing the top-N co-occurring
    domain pairs as a colour-filled bar chart.

    The file can be rendered with:
        xmgrace domain_cooccurrence.agr
    or non-interactively:
        gracebat domain_cooccurrence.agr -hdevice PNG \\
                 -hardcopy -printfile output.png

    Parameters
    ----------
    cooccurrence : Counter from compute_domain_cooccurrence()
    output_path  : destination .agr file
    top_n        : number of domain pairs to include
    title        : graph title string

    Returns
    -------
    str  — full .agr file content
    """
    # Off-diagonal pairs only (true co-occurrence, not frequency)
    pairs = [(p, c) for p, c in cooccurrence.items() if p[0] != p[1]]
    pairs.sort(key=lambda x: x[1], reverse=True)
    top = pairs[:top_n]

    if not top:
        # fall back to single-domain frequency when no pairs exist
        top = [
            (p, c) for p, c in cooccurrence.most_common(top_n)
            if p[0] == p[1]
        ]

    n         = len(top)
    max_count = max(c for _, c in top) if top else 10
    y_max     = max_count * 1.18

    # Shortened labels for axis ticks
    labels = []
    for (d1, d2), _ in top:
        if d1 == d2:
            lbl = d1[:22]
        else:
            lbl = f"{d1[:14]}/{d2[:14]}"
        labels.append(lbl)

    # ── header ──────────────────────────────────────────────────────────────
    lines = [
        "# Grace project file",
        "# Generated by Protein Domain Analysis Toolkit",
        f"# Title: {title}",
        "#",
        "@version 50125",
        "@page size 1200, 800",
        "@default linewidth 2.0",
        "@default charsize 1.200000",
        "@background color 0",
        "@page background fill on",
        "@timestamp off",
        "",
        _AGR_FONT_MAP.strip(),
        "",
        _AGR_COLOR_MAP.strip(),
        "",
        "# ── graph ──────────────────────────────────────────────────────",
        "@g0 on",
        "@g0 hidden false",
        "@g0 type XY",
        "@with g0",
        f'@    world -0.5, 0, {n + 0.5:.1f}, {y_max:.1f}',
        "@    view 0.150000, 0.200000, 1.100000, 0.900000",
        f'@    title "{title}"',
        "@    title font 2",
        "@    title size 1.750000",
        "@    title color 1",
        '@    subtitle "Source: UniProt (Swiss-Prot reviewed)"',
        "@    subtitle font 0",
        "@    subtitle size 1.000000",
        "@    subtitle color 7",
        "",
        "# ── x-axis ──────────────────────────────────────────────────────",
        "@    xaxis on",
        "@    xaxis bar on",
        "@    xaxis bar color 1",
        "@    xaxis bar linewidth 2.0",
        '@    xaxis label "Domain Pair (ranked by frequency)"',
        "@    xaxis label char size 1.300000",
        "@    xaxis label font 0",
        "@    xaxis label color 1",
        "@    xaxis tick on",
        f"@    xaxis tick major {max(1, n // 5)}",
        "@    xaxis tick minor ticks 0",
        "@    xaxis ticklabel on",
        "@    xaxis ticklabel type spec",
        "@    xaxis tick spec type both",
        f"@    xaxis tick spec {n}",
    ]

    for i, lbl in enumerate(labels):
        lines.append(f'@    xaxis tick major {i}, {i + 1}')
        lines.append(f'@    xaxis ticklabel {i}, "{lbl}"')

    lines += [
        "@    xaxis ticklabel char size 0.750000",
        "@    xaxis ticklabel angle 45",
        "@    xaxis ticklabel font 0",
        "",
        "# ── y-axis ──────────────────────────────────────────────────────",
        "@    yaxis on",
        "@    yaxis bar on",
        "@    yaxis bar color 1",
        "@    yaxis bar linewidth 2.0",
        '@    yaxis label "Co-occurrence Count"',
        "@    yaxis label char size 1.300000",
        "@    yaxis label font 0",
        "@    yaxis label color 1",
        "@    yaxis tick on",
        f"@    yaxis tick major {max(1, int(max_count / 6))}",
        "@    yaxis tick minor ticks 4",
        "@    yaxis ticklabel char size 1.000000",
        "",
        "# ── dataset ─────────────────────────────────────────────────────",
        "@    s0 type bar",
        "@    s0 color 4",
        "@    s0 fill color 4",
        "@    s0 fill pattern 1",
        "@    s0 line color 1",
        "@    s0 line linewidth 1.0",
        "@    s0 basetype zero",
        "@    s0 avalue on",
        "@    s0 avalue type 2",
        "@    s0 avalue char size 0.800000",
        "@    s0 avalue color 1",
        "@    s0 avalue font 0",
        "@    s0 avalue rot 0",
        "@    s0 avalue format decimal",
        "@    s0 avalue prec 0",
        '@    s0 legend "Co-occurrence frequency"',
        "@    legend on",
        "@    legend loctype view",
        "@    legend 0.75, 0.92",
        "@    legend char size 1.000000",
        "@    frame type 0",
        "@target G0.S0",
        "@type bar",
    ]

    for i, (_, cnt) in enumerate(top):
        lines.append(f"{i + 1} {cnt}")

    lines += [
        "&",
        "",
        "# ── domain pair annotation strings ─────────────────────────────",
    ]

    for i, ((d1, d2), cnt) in enumerate(top):
        lbl = f"{d1} + {d2}" if d1 != d2 else d1
        lines.append(f'# pair {i + 1}: {lbl} (n={cnt})')

    content = "\n".join(lines) + "\n"

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(content)
    print(f"  [XMGrace] .agr file → {out}")
    return content


# ════════════════════════════════════════════════════════════════════════════
#  4. BED → FASTA  (shell: Convert BED coordinate to FASTA.sh)
# ════════════════════════════════════════════════════════════════════════════

def bed_to_fasta(bed_file: str, genome_dir: str, output_fasta: str) -> int:
    """
    Extract FASTA sequences for BED coordinate ranges.

    Uses pyfaidx when available; falls back to invoking the legacy
    `fastacmd` binary from NCBI toolkit (original script behaviour).

    Parameters
    ----------
    bed_file     : path to BED3+ file
    genome_dir   : directory containing per-chromosome FASTA files
                   named <chrom>.fa
    output_fasta : path for the output FASTA file

    Returns
    -------
    int  number of records written
    """
    records_written = 0
    entries: list[tuple[str, int, int]] = []

    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            entries.append((parts[0], int(parts[1]), int(parts[2])))

    with open(output_fasta, "w") as out:
        for chrom, start, end in entries:
            fa_file = os.path.join(genome_dir, f"{chrom}.fa")

            if HAS_PYFAIDX and os.path.exists(fa_file):
                genome = Fasta(fa_file)
                if chrom in genome:
                    seq = genome[chrom][start:end].seq
                    out.write(f">{chrom}:{start}-{end}\n{seq}\n")
                    records_written += 1
            else:
                # original fastacmd fallback
                cmd = [
                    "fastacmd", "-d", fa_file, "-p", "F",
                    "-L", f"{start},{end}", "-s", chrom, "-S", "1",
                ]
                try:
                    result = subprocess.run(
                        cmd, capture_output=True, text=True, timeout=30
                    )
                    if result.stdout:
                        out.write(result.stdout)
                        records_written += 1
                except (subprocess.SubprocessError, FileNotFoundError):
                    print(f"  [BED→FASTA] Warning: could not extract "
                          f"{chrom}:{start}-{end}")

    print(f"  [BED→FASTA] {records_written} records → {output_fasta}")
    return records_written


# ════════════════════════════════════════════════════════════════════════════
#  5. Remove redundant FASTA  (Perl: Remove Redundant FASTA.pl)
# ════════════════════════════════════════════════════════════════════════════

def remove_redundant_fasta(input_fasta: str, output_fasta: str) -> tuple[int, int]:
    """
    Collapse sequences sharing identical residue strings.

    When multiple proteins share the same sequence, their IDs are
    concatenated with ';', matching the Perl script's behaviour.

    Requires BioPython.

    Returns
    -------
    (total_in, unique_out)
    """
    if not HAS_BIOPYTHON:
        raise RuntimeError("BioPython is required.  pip install biopython")

    seen:    dict[str, SeqRecord] = {}
    headers: dict[str, list[str]] = {}
    total = 0

    with open(input_fasta) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            total += 1
            seq_str = str(rec.seq).upper()
            if seq_str in seen:
                headers[seq_str].append(rec.id)
            else:
                seen[seq_str]    = rec
                headers[seq_str] = [rec.id]

    unique = 0
    with open(output_fasta, "w") as out:
        for seq_str, rec in seen.items():
            ids = headers[seq_str]
            if len(ids) > 1:
                rec.id          = ";".join(ids)
                rec.description = f"redundant_group n={len(ids)}"
            SeqIO.write(rec, out, "fasta")
            unique += 1

    print(f"  [RedFASTA] {total} in → {unique} unique out "
          f"({total - unique} removed)")
    return total, unique


# ════════════════════════════════════════════════════════════════════════════
#  6. Biological network duplicate detection
#     (Duplicates in Biological Network.sh)
# ════════════════════════════════════════════════════════════════════════════

def detect_network_duplicates(
    input_file:  str,
    output_file: str,
    col1:        int = 0,
    col2:        int = 1,
    col3:        int = 5,
    sep:         str = "\t",
) -> pd.DataFrame:
    """
    Count and remove duplicate edges in a biological network file.

    Edge (A, interaction, B) is treated as identical to (B, interaction, A).
    The canonical form always places the lexicographically smaller node first.

    The shell script used an awk pipeline:
        awk '$1>$3{next}{print $2,$1,$3}' | sort | uniq -c

    Parameters
    ----------
    col1 / col3 : column indices for the two interacting nodes
    col2        : column index for the interaction type / score
    """
    df = pd.read_csv(input_file, sep=sep, header=None, dtype=str)

    nc = len(df.columns)
    c1 = min(col1, nc - 1)
    c2 = min(col2, nc - 1)
    c3 = min(col3, nc - 1)

    def _norm(row: pd.Series) -> tuple[str, str, str]:
        a, b = str(row.iloc[c1]), str(row.iloc[c3])
        if a > b:
            a, b = b, a
        return (a, str(row.iloc[c2]), b)

    df["_key"] = df.apply(_norm, axis=1)

    counts = df["_key"].value_counts()
    result = (
        df.drop_duplicates(subset=["_key"])
          .copy()
          .assign(duplicate_count=lambda x: x["_key"].map(counts))
          .drop(columns=["_key"])
    )
    result.to_csv(output_file, sep=sep, index=False, header=False)

    dups = (counts > 1).sum()
    print(f"  [DupNet] {dups} duplicate edges found in "
          f"{len(counts)} total pairs")
    return result


# ════════════════════════════════════════════════════════════════════════════
#  7. Overlapping domain region collapse
#     (Domain cooccurence.sh / Domain cooccurence-o.sh core algorithm)
# ════════════════════════════════════════════════════════════════════════════

def process_overlap_regions(
    input_file:    str,
    output_min:    str = "min_nmbr.txt",
    output_max:    str = "max_nmbr.txt",
    output_region: str = "region_limit.txt",
) -> pd.DataFrame:
    """
    Collapse overlapping BLAST/domain hit regions per protein.

    Algorithm (from Domain cooccurence-o.sh):
    - Sort hits by (protein, begin)
    - For each protein, merge overlapping intervals
    - Within a cluster, retain the hit with lowest e-value (best hit)
    - Emit min-begin and max-end boundary tables

    Expected input columns (0-based):
        0  protein_id
        3  begin
        4  end
        5  e-value exponent

    Returns
    -------
    pd.DataFrame of non-overlapping best hits
    """
    try:
        df = pd.read_csv(input_file, sep="\t", header=None, dtype=str)
    except Exception:
        df = pd.read_csv(input_file, sep=r"\s+", header=None, dtype=str)

    nc = len(df.columns)
    c_prot = 0
    c_beg  = min(3, nc - 1)
    c_end  = min(4, nc - 1)
    c_eval = min(5, nc - 1)

    df[c_beg]  = pd.to_numeric(df[c_beg],  errors="coerce").fillna(0).astype(int)
    df[c_end]  = pd.to_numeric(df[c_end],  errors="coerce").fillna(0).astype(int)
    df[c_eval] = pd.to_numeric(df[c_eval], errors="coerce").fillna(999.0)

    # min begin / max end per protein  (awk one-liners from original script)
    grp = df.groupby(df[c_prot])
    grp[c_beg].min().to_csv(output_min, sep=" ", header=False)
    grp[c_end].max().to_csv(output_max, sep=" ", header=False)

    region_df = pd.concat(
        [grp[c_beg].min(), grp[c_end].max()], axis=1
    )
    region_df.columns = ["min_begin", "max_end"]
    region_df.to_csv(output_region, sep="\t")

    # cluster + best-hit selection
    df_sorted = df.sort_values([c_prot, c_beg]).reset_index(drop=True)

    clusters: list[pd.Series] = []
    cur_prot: str | None = None
    cur_end  = 0
    cur_best_eval = float("inf")
    cur_best: pd.Series | None = None

    for _, row in df_sorted.iterrows():
        prot  = str(row[c_prot])
        begin = int(row[c_beg])
        end   = int(row[c_end])
        eval_ = float(row[c_eval])

        if prot == cur_prot and begin <= cur_end:
            # extend cluster
            cur_end = max(cur_end, end)
            if eval_ < cur_best_eval:
                cur_best_eval = eval_
                cur_best      = row
        else:
            if cur_best is not None:
                clusters.append(cur_best)
            cur_prot      = prot
            cur_end       = end
            cur_best_eval = eval_
            cur_best      = row

    if cur_best is not None:
        clusters.append(cur_best)

    result = pd.DataFrame(clusters).reset_index(drop=True)
    print(f"  [OverlapCollapse] {len(df)} rows → "
          f"{len(result)} non-overlapping clusters")
    return result


# ════════════════════════════════════════════════════════════════════════════
#  8. Visualisation helper
# ════════════════════════════════════════════════════════════════════════════

def plot_cooccurrence_heatmap(
    matrix:     pd.DataFrame,
    output_path: Path | str,
    title:       str = "Domain Co-occurrence Heatmap",
    top_n:       int = 25,
) -> None:
    """
    Render domain co-occurrence as a seaborn heatmap (PNG, 150 dpi).

    Only the top_n most frequent domains are shown to keep the plot legible.
    """
    if matrix.empty:
        print("  [Heatmap] Empty matrix — skipping")
        return

    freq   = matrix.sum(axis=1).sort_values(ascending=False)
    top    = freq.head(top_n).index.tolist()
    sub    = matrix.loc[top, top]

    fig, ax = plt.subplots(figsize=(14, 12))
    sns.heatmap(
        sub,
        ax=ax,
        cmap="YlOrRd",
        annot=(len(top) <= 15),
        fmt="d",
        linewidths=0.4,
        cbar_kws={"label": "Co-occurrence count", "shrink": 0.8},
    )
    ax.set_title(title, fontsize=14, fontweight="bold", pad=16)
    ax.set_xlabel("Domain", fontsize=11)
    ax.set_ylabel("Domain", fontsize=11)
    plt.xticks(rotation=45, ha="right", fontsize=7)
    plt.yticks(rotation=0,  fontsize=7)
    plt.tight_layout()

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  [Heatmap] → {out}")


# ════════════════════════════════════════════════════════════════════════════
#  9. Interactive CLI menu
# ════════════════════════════════════════════════════════════════════════════

MENU = {
    "1": "Mine UniProt proteomics data",
    "2": "Compute domain co-occurrence from TSV",
    "3": "Generate XMGrace .agr output",
    "4": "Convert BED coordinates to FASTA",
    "5": "Remove redundant FASTA sequences",
    "6": "Detect duplicates in biological network",
    "7": "Collapse overlapping domain regions",
    "8": "Run full S. cerevisiae workflow  [recommended demo]",
    "9": "Exit",
}


def _banner():
    print("""
╔══════════════════════════════════════════════════════════╗
║       Protein Domain Analysis Toolkit  v2.0             ║
║       Dr. Arli Aditya Parikesit  |  i3L University      ║
╚══════════════════════════════════════════════════════════╝""")


def _print_menu():
    _banner()
    for k, v in MENU.items():
        print(f"  [{k}] {v}")
    print()


def _ask(prompt: str, default: str = "") -> str:
    val = input(f"  {prompt} [{default}]: ").strip()
    return val if val else default


def interactive_menu():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    while True:
        _print_menu()
        choice = input("  Select option: ").strip()

        # ── 1: Mine UniProt ────────────────────────────────────────────────
        if choice == "1":
            org_id   = int(_ask("Organism Tax ID",     str(SCEREVISIAE_TAXID)))
            max_r    = int(_ask("Max proteins",         "300"))
            reviewed = _ask("Swiss-Prot only? (y/n)", "y").lower() == "y"
            mine_uniprot(org_id, reviewed=reviewed, max_results=max_r)

        # ── 2: Co-occurrence ───────────────────────────────────────────────
        elif choice == "2":
            tsv = _ask("Input TSV path")
            if not Path(tsv).exists():
                print(f"  File not found: {tsv}"); continue
            df = pd.read_csv(tsv, sep="\t")
            matrix, cooc = compute_domain_cooccurrence(df)
            if matrix.empty:
                print("  No domain data found."); continue
            print(f"\n  Matrix: {matrix.shape[0]} × {matrix.shape[1]} domains")
            print("  Top-10 co-occurring pairs:")
            for pair, cnt in [
                (p, c) for p, c in cooc.most_common(20) if p[0] != p[1]
            ][:10]:
                print(f"    {pair[0]}  +  {pair[1]}  →  {cnt}")
            out = RESULTS_DIR / "cooccurrence_matrix.tsv"
            matrix.to_csv(out, sep="\t")
            print(f"  Matrix saved → {out}")

        # ── 3: XMGrace ─────────────────────────────────────────────────────
        elif choice == "3":
            mat_path = RESULTS_DIR / "cooccurrence_matrix.tsv"
            if not mat_path.exists():
                print("  Run option 2 first."); continue
            mat = pd.read_csv(mat_path, sep="\t", index_col=0)
            cooc: Counter = Counter()
            for d1 in mat.index:
                for d2 in mat.columns:
                    cnt = int(mat.loc[d1, d2])
                    if cnt > 0 and str(d1) <= str(d2):
                        cooc[(str(d1), str(d2))] = cnt
            top_n = int(_ask("Top N pairs in plot", "20"))
            generate_xmgrace_cooccurrence(
                cooc, RESULTS_DIR / "domain_cooccurrence.agr", top_n=top_n
            )

        # ── 4: BED → FASTA ─────────────────────────────────────────────────
        elif choice == "4":
            bed  = _ask("BED file path")
            gdir = _ask("Genome directory")
            out  = _ask("Output FASTA path", "output.fasta")
            if Path(bed).exists():
                bed_to_fasta(bed, gdir, out)
            else:
                print(f"  File not found: {bed}")

        # ── 5: Remove redundant FASTA ──────────────────────────────────────
        elif choice == "5":
            in_fa  = _ask("Input FASTA path")
            out_fa = _ask("Output FASTA path", "nonredundant.fa")
            if Path(in_fa).exists():
                remove_redundant_fasta(in_fa, out_fa)
            else:
                print(f"  File not found: {in_fa}")

        # ── 6: Network duplicates ─────────────────────────────────────────
        elif choice == "6":
            in_f  = _ask("Input network file")
            out_f = _ask("Output file", "deduped_network.tsv")
            if Path(in_f).exists():
                detect_network_duplicates(in_f, out_f)
            else:
                print(f"  File not found: {in_f}")

        # ── 7: Overlap collapse ───────────────────────────────────────────
        elif choice == "7":
            in_f = _ask("Input overlap list file")
            if Path(in_f).exists():
                process_overlap_regions(in_f)
            else:
                print(f"  File not found: {in_f}")

        # ── 8: Full S. cerevisiae workflow ────────────────────────────────
        elif choice == "8":
            _run_scerevisiae_workflow()

        # ── 9: Exit ───────────────────────────────────────────────────────
        elif choice == "9":
            print("\n  Farewell.\n")
            sys.exit(0)

        else:
            print(f"  Unknown option: {choice}")


# ════════════════════════════════════════════════════════════════════════════
#  10. Full S. cerevisiae demonstration workflow
# ════════════════════════════════════════════════════════════════════════════

def _run_scerevisiae_workflow(max_proteins: int = 300, top_n: int = 20) -> dict:
    """
    End-to-end pipeline:
      mine UniProt → parse domains → co-occurrence matrix →
      XMGrace .agr bar chart → PNG heatmap

    Use-case: S. cerevisiae (Baker's yeast) reviewed proteome.
    Yeast is ideal as a benchmark because its ~6,000 proteins are
    exhaustively curated in Swiss-Prot and span all major eukaryotic
    domain families.

    Returns
    -------
    dict with keys: df, matrix, cooccurrence, agr_path, heatmap_path
    """
    _banner()
    print("\n  ── S. cerevisiae Domain Co-occurrence Workflow ──\n")

    # 1. Mine
    print("  [Step 1/4] Mining UniProt ...")
    df = mine_uniprot(
        SCEREVISIAE_TAXID, reviewed=True, max_results=max_proteins
    )

    # 2. Co-occurrence
    print("\n  [Step 2/4] Computing domain co-occurrence ...")
    matrix, cooc = compute_domain_cooccurrence(df)

    if matrix.empty:
        print("  No domain annotations in retrieved data.  "
              "Try increasing max_proteins or check UniProt connectivity.")
        return {"df": df, "matrix": matrix, "cooccurrence": cooc}

    mat_path = RESULTS_DIR / "scerevisiae_cooccurrence_matrix.tsv"
    matrix.to_csv(mat_path, sep="\t")
    print(f"  Matrix ({matrix.shape[0]} domains) → {mat_path}")

    print("\n  Top-15 co-occurring domain pairs:")
    for pair, cnt in [
        (p, c) for p, c in cooc.most_common(25) if p[0] != p[1]
    ][:15]:
        print(f"    {pair[0]:<35}  +  {pair[1]:<35}  n={cnt}")

    # 3. XMGrace
    print("\n  [Step 3/4] Generating XMGrace output ...")
    agr_path = RESULTS_DIR / "scerevisiae_domain_cooccurrence.agr"
    generate_xmgrace_cooccurrence(
        cooc,
        agr_path,
        top_n=top_n,
        title="Domain Co-occurrence in Saccharomyces cerevisiae",
    )
    print(f"  Render with:  xmgrace {agr_path}")
    print(f"  Or PNG batch: gracebat {agr_path} "
          "-hdevice PNG -hardcopy -printfile scerevisiae_cooc.png")

    # 4. Heatmap
    print("\n  [Step 4/4] Generating co-occurrence heatmap ...")
    hm_path = RESULTS_DIR / "scerevisiae_domain_heatmap.png"
    plot_cooccurrence_heatmap(
        matrix, hm_path,
        title="Domain Co-occurrence — Saccharomyces cerevisiae "
              "(UniProt Swiss-Prot reviewed)",
        top_n=25,
    )

    print("\n  ── Workflow complete ──")
    print(f"  Results in: {RESULTS_DIR.resolve()}/\n")

    return {
        "df":           df,
        "matrix":       matrix,
        "cooccurrence": cooc,
        "agr_path":     agr_path,
        "heatmap_path": hm_path,
    }


# ════════════════════════════════════════════════════════════════════════════
#  CLI entry-point
# ════════════════════════════════════════════════════════════════════════════

def main():
    import argparse

    parser = argparse.ArgumentParser(
        prog="protein_domain_toolkit",
        description="Protein Domain Analysis Toolkit — Dr. Arli Aditya Parikesit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  python protein_domain_toolkit.py               # interactive menu
  python protein_domain_toolkit.py --workflow    # S. cerevisiae demo
  python protein_domain_toolkit.py --mine --org 559292 --max 500
  python protein_domain_toolkit.py --cooccurrence results/uniprot_559292.tsv
  python protein_domain_toolkit.py --xmgrace results/cooccurrence_matrix.tsv
""",
    )
    parser.add_argument("--workflow",       action="store_true",
                        help="Run full S. cerevisiae workflow")
    parser.add_argument("--mine",           action="store_true",
                        help="Mine UniProt only")
    parser.add_argument("--org",            type=int, default=SCEREVISIAE_TAXID,
                        metavar="TAXID",
                        help="Organism taxonomy ID (default 559292)")
    parser.add_argument("--max",            type=int, default=300,
                        metavar="N",
                        help="Max proteins to retrieve (default 300)")
    parser.add_argument("--cooccurrence",   metavar="TSV",
                        help="Compute co-occurrence from TSV file")
    parser.add_argument("--xmgrace",        metavar="MATRIX_TSV",
                        help="Generate XMGrace .agr from co-occurrence matrix")
    parser.add_argument("--bed-to-fasta",   nargs=3,
                        metavar=("BED", "GENOME_DIR", "OUT_FA"),
                        help="Convert BED to FASTA")
    parser.add_argument("--remove-redundant", nargs=2,
                        metavar=("IN_FA", "OUT_FA"),
                        help="Remove redundant FASTA sequences")
    parser.add_argument("--detect-duplicates", nargs=2,
                        metavar=("IN", "OUT"),
                        help="Detect duplicates in biological network")

    args = parser.parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # No flags → interactive
    if len(sys.argv) == 1:
        interactive_menu()
        return

    if args.workflow:
        _run_scerevisiae_workflow(max_proteins=args.max)
        return

    if args.mine:
        mine_uniprot(args.org, max_results=args.max)

    if args.cooccurrence:
        df = pd.read_csv(args.cooccurrence, sep="\t")
        matrix, cooc = compute_domain_cooccurrence(df)
        matrix.to_csv(RESULTS_DIR / "cooccurrence_matrix.tsv", sep="\t")
        generate_xmgrace_cooccurrence(
            cooc, RESULTS_DIR / "domain_cooccurrence.agr"
        )
        plot_cooccurrence_heatmap(
            matrix, RESULTS_DIR / "domain_heatmap.png"
        )

    if args.xmgrace:
        mat = pd.read_csv(args.xmgrace, sep="\t", index_col=0)
        cooc: Counter = Counter()
        for d1 in mat.index:
            for d2 in mat.columns:
                cnt = int(mat.loc[d1, d2])
                if cnt > 0 and str(d1) <= str(d2):
                    cooc[(str(d1), str(d2))] = cnt
        generate_xmgrace_cooccurrence(
            cooc, RESULTS_DIR / "domain_cooccurrence.agr"
        )

    if args.bed_to_fasta:
        bed_to_fasta(*args.bed_to_fasta)

    if args.remove_redundant:
        remove_redundant_fasta(*args.remove_redundant)

    if args.detect_duplicates:
        detect_network_duplicates(*args.detect_duplicates)


if __name__ == "__main__":
    main()
