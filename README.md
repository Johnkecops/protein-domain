# Protein Domain Analysis Toolkit

> **Integrating legacy bioinformatics shell scripts into a reproducible Python pipeline with a Streamlit web interface.**

**Author:** Dr.rer.nat. Arli Aditya Parikesit, S.Si., M.Si.  
**Affiliation:** Department of Biotechnology, School of Health and Life Sciences, i3L University, Jakarta  
**ORCID:** [0000-0001-8716-3926](https://orcid.org/0000-0001-8716-3926)  
**PhD Dissertation:** [DOI 10.6084/M9.FIGSHARE.964089](https://doi.org/10.6084/M9.FIGSHARE.964089)

---

## Version History

| Version | Date       | Changes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| ------- | ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **1.1** | 2026-05-10 | Code review fixes: 21 issues resolved across `protein_domain_toolkit.py` and `streamlit_app.py`. Medium severity: silent domain-drop logging, explicit `ft_domain` column detection, XMGrace label disambiguation and fallback warning, diagonal-corrected heatmap ranking, complete return dict from workflow, `--mine` CLI tip, genome dir validation in Streamlit, stale session state warning. Low severity: `UNIPROT_PAGE_SIZE` constant, 3-attempt retry on UniProt requests, narrowed coordinate filter in `_parse_domains`, empty DataFrame contract comment, `@version` annotation, annotation threshold comment, `--workflow` flag conflict warning, deprecated `use_column_width` replaced, `st.empty()` log overwrite fixed, placeholder image removed. See `recap.md` for full issue catalogue. |
| **1.0** | 2026-05-01 | Initial release. Consolidated five legacy shell/Perl scripts into unified Python toolkit and Streamlit web app. Live S. cerevisiae use-case via UniProt REST API. XMGrace `.agr` export and seaborn heatmap output.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

---

## Background

This repository archives and modernises shell and Perl scripts originally developed during a PhD project (2008–2012) on evolutionary protein domain distributions across all three domains of life.  The scripts have been consolidated into a single Python toolkit and a Streamlit web application, with a live *Saccharomyces cerevisiae* use-case demonstrating the full pipeline.

### Original Publications

| Year | Reference |
|------|-----------|
| 2014 | Parikesit, A.A., Steiner, L., Stadler, P.F., Prohaska, S.J. — *Malaysian J. Fund. Appl. Sci.*, 10(2), 65–75. [DOI](https://doi.org/10.11113/mjfas.v10n2.57) |
| 2011 | Parikesit, A.A., Stadler, P.F., Prohaska, S.J. — *Genes*, 2(4), 912–924. [DOI](https://doi.org/10.3390/genes2040912) |
| 2010 | Parikesit, A.A., Stadler, P.F., Prohaska, S.J. — *Ger. Conf. Bioinform.*, P-173, 93–102. [DOI](https://doi.org/10.13140/2.1.2036.2248) |

---

## Repository Structure

```
protein-domain-main/
├── protein_domain_toolkit.py     # Integrated CLI toolkit (menu-driven)
├── streamlit_app.py              # Streamlit web application
├── requirements.txt              # Python dependencies
├── README.md                     # This file
├── LICENSE.txt                   # MIT licence
├── results/                      # Auto-created output directory
│   ├── uniprot_559292.tsv        # Mined S. cerevisiae proteome
│   ├── scerevisiae_cooccurrence_matrix.tsv
│   ├── scerevisiae_domain_cooccurrence.agr   # XMGrace project file
│   └── scerevisiae_domain_heatmap.png
│
│   [Legacy scripts — preserved for reference]
├── Domain cooccurence.sh
├── Domain cooccurence-o.sh
├── Convert BED coordinate to FASTA.sh
├── Duplicates in Biological Network.sh
└── Remove Redundant FASTA.pl
```

---

## Tools Integrated

| # | Legacy Script | Python Function | Description |
|---|---------------|-----------------|-------------|
| 1 | `Domain cooccurence.sh` | `process_overlap_regions()` | Collapse overlapping domain hit clusters; retain best e-value hit |
| 2 | `Domain cooccurence-o.sh` | `process_overlap_regions()` | Extended overlap logic with min/max boundary tables |
| 3 | `Convert BED coordinate to FASTA.sh` | `bed_to_fasta()` | Extract FASTA sequences from BED genomic coordinates |
| 4 | `Duplicates in Biological Network.sh` | `detect_network_duplicates()` | Normalise and count duplicate edges in network files |
| 5 | `Remove Redundant FASTA.pl` | `remove_redundant_fasta()` | Collapse identical sequences; concatenate headers |
| 6 | *(new)* | `mine_uniprot()` | Mine UniProt REST API for any organism's reviewed proteome |
| 7 | *(new)* | `compute_domain_cooccurrence()` | Build pairwise domain co-occurrence matrix |
| 8 | *(new)* | `generate_xmgrace_cooccurrence()` | Write XMGrace `.agr` bar-chart project file |
| 9 | *(new)* | `plot_cooccurrence_heatmap()` | Render seaborn heatmap (PNG) |

---

## Installation

```bash
# Clone
git clone https://github.com/arliadityaparikesit/protein-domain.git
cd protein-domain-main

# Create environment (conda recommended)
conda create -n proteindomain python=3.11
conda activate proteindomain

# Install dependencies
pip install -r requirements.txt
```

---

## Quick Start

### Interactive CLI (recommended for command-line users)

```bash
python protein_domain_toolkit.py
```

The menu presents all eight tools:

```
╔══════════════════════════════════════════════════════════╗
║       Protein Domain Analysis Toolkit  v2.0             ║
║       Dr. Arli Aditya Parikesit  |  i3L University      ║
╚══════════════════════════════════════════════════════════╝
  [1] Mine UniProt proteomics data
  [2] Compute domain co-occurrence from TSV
  [3] Generate XMGrace .agr output
  [4] Convert BED coordinates to FASTA
  [5] Remove redundant FASTA sequences
  [6] Detect duplicates in biological network
  [7] Collapse overlapping domain regions
  [8] Run full S. cerevisiae workflow  [recommended demo]
  [9] Exit
```

### One-command S. cerevisiae workflow

```bash
python protein_domain_toolkit.py --workflow
```

This runs all four pipeline stages automatically (mine → analyse → XMGrace → heatmap) and writes results to `results/`.

### Streamlit web app

```bash
streamlit run streamlit_app.py
```

Open [http://localhost:8501](http://localhost:8501) in your browser.

---

## S. cerevisiae Use-Case

*Saccharomyces cerevisiae* (Baker's yeast, TaxID 559292) is used as the canonical demonstration organism because:

- Its proteome (~6,000 proteins) is exhaustively curated in UniProt Swiss-Prot
- It spans all major eukaryotic domain families documented in InterPro/Pfam
- It served as a primary benchmark organism in the original PhD research
- Results are immediately reproducible by any researcher with internet access

### What the workflow produces

1. **`uniprot_559292.tsv`** — Reviewed proteome with `ft_domain` annotations retrieved via the UniProt REST API
2. **`scerevisiae_cooccurrence_matrix.tsv`** — Symmetric N×N matrix; cell [i,j] = number of proteins carrying both domain i and domain j
3. **`scerevisiae_domain_cooccurrence.agr`** — XMGrace project file; bar chart of the top-20 co-occurring domain pairs
4. **`scerevisiae_domain_heatmap.png`** — Seaborn heatmap of the top-25 most frequent domains

### XMGrace rendering

```bash
# Interactive
xmgrace results/scerevisiae_domain_cooccurrence.agr

# Non-interactive PNG export (requires gracebat)
gracebat results/scerevisiae_domain_cooccurrence.agr \
  -hdevice PNG -hardcopy -printfile scerevisiae_cooc.png

# EPS for LaTeX
gracebat results/scerevisiae_domain_cooccurrence.agr \
  -hdevice EPS -hardcopy -printfile scerevisiae_cooc.eps
```

---

## CLI Reference

```
usage: protein_domain_toolkit [-h] [--workflow] [--mine] [--org TAXID]
                               [--max N] [--cooccurrence TSV]
                               [--xmgrace MATRIX_TSV]
                               [--bed-to-fasta BED GENOME_DIR OUT_FA]
                               [--remove-redundant IN_FA OUT_FA]
                               [--detect-duplicates IN OUT]

Options
  --workflow            Run full S. cerevisiae workflow
  --mine                Mine UniProt data only
  --org TAXID           Organism taxonomy ID (default 559292)
  --max N               Max proteins to retrieve (default 300)
  --cooccurrence TSV    Compute co-occurrence from TSV
  --xmgrace MATRIX_TSV  Generate XMGrace .agr from matrix
  --bed-to-fasta        Convert BED → FASTA
  --remove-redundant    Remove redundant FASTA sequences
  --detect-duplicates   Detect duplicate edges in network
```

---

## Algorithm Notes

### Domain co-occurrence (core algorithm from PhD research)

For each reviewed protein with ≥2 annotated domains:
1. Extract unique domain names from the `ft_domain` field
2. Enumerate all unordered pairs (d_i, d_j) where i < j
3. Increment the co-occurrence counter for each pair
4. Build a symmetric matrix from the counter

The matrix entry M[i,j] represents how many proteins in the proteome carry both domain i and domain j simultaneously — a direct measure of functional coupling and multi-domain architecture.

### Overlap collapse (from Domain cooccurence*.sh)

Given a sorted list of BLAST domain hits per protein:
1. Sort by (protein_id, begin_position)
2. Merge intervals with begin ≤ current_end (overlapping)
3. Within each cluster, retain the hit with the lowest e-value
4. Output min-begin and max-end boundary tables

This eliminates redundant overlapping hits before co-occurrence counting.

---

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| biopython | ≥1.81 | FASTA I/O, sequence handling |
| pandas | ≥2.0 | Data frames, TSV I/O |
| numpy | ≥1.24 | Numerical operations |
| matplotlib | ≥3.7 | Base plotting |
| seaborn | ≥0.12 | Heatmap rendering |
| requests | ≥2.31 | UniProt REST API |
| pyfaidx | ≥0.7 | Indexed FASTA access (BED→FASTA) |
| streamlit | ≥1.32 | Web interface |

---

## Licence

MIT — see [LICENSE.txt](LICENSE.txt).

---

## Citation

If this toolkit contributes to published work, please cite:

```bibtex
@article{parikesit2014pitfalls,
  author  = {Parikesit, Arli Aditya and Steiner, Lydia and
             Stadler, Peter F. and Prohaska, Sonja J.},
  title   = {Pitfalls of Ascertainment Biases in Genome Annotations —
             Computing Comparable Protein Domain Distributions in Eukarya},
  journal = {Malaysian Journal of Fundamental and Applied Sciences},
  volume  = {10},
  number  = {2},
  pages   = {65--75},
  year    = {2014},
  doi     = {10.11113/mjfas.v10n2.57}
}

@article{parikesit2011evolution,
  author  = {Parikesit, Arli Aditya and Stadler, Peter F. and
             Prohaska, Sonja J.},
  title   = {Evolution and Quantitative Comparison of Genome-Wide
             Protein Domain Distributions},
  journal = {Genes},
  volume  = {2},
  number  = {4},
  pages   = {912--924},
  year    = {2011},
  doi     = {10.3390/genes2040912}
}
```

**AI Assistance Disclaimer**: This codebase was developed with the assistance of Claude Code. While the AI provided code generation, debugging, and structural support, the human developer maintains full responsibility for reviewing, testing, and maintaining all content and functionality.
