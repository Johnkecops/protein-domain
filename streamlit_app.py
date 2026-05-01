#!/usr/bin/env python3
"""
Protein Domain Analysis Toolkit — Streamlit Web App
====================================================
Author  : Dr. Arli Aditya Parikesit
Affil.  : Department of Bioinformatics, i3L University, Jakarta
ORCID   : 0000-0001-8716-3926

Run locally:
    streamlit run streamlit_app.py

The app mirrors every feature in protein_domain_toolkit.py through
an interactive browser UI — no command-line knowledge required.
"""

from __future__ import annotations

import io
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

# Ensure toolkit module resolves when app is started from project root
sys.path.insert(0, str(Path(__file__).parent))

from protein_domain_toolkit import (
    RESULTS_DIR,
    SCEREVISIAE_TAXID,
    _parse_domains,
    bed_to_fasta,
    compute_domain_cooccurrence,
    detect_network_duplicates,
    generate_xmgrace_cooccurrence,
    mine_uniprot,
    plot_cooccurrence_heatmap,
    process_overlap_regions,
    remove_redundant_fasta,
)

# ── page configuration ───────────────────────────────────────────────────────
st.set_page_config(
    page_title="Protein Domain Analysis Toolkit",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "Get help": "https://orcid.org/0000-0001-8716-3926",
        "About": (
            "Protein Domain Analysis Toolkit\n"
            "Dr. Arli Aditya Parikesit | i3L University Jakarta\n"
            "ORCID: 0000-0001-8716-3926"
        ),
    },
)

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── custom CSS ────────────────────────────────────────────────────────────────
st.markdown(
    """
<style>
    .block-container { padding-top: 1.5rem; }
    .sidebar-content { font-size: 0.9rem; }
    h1 { font-size: 1.6rem !important; }
    h2 { font-size: 1.25rem !important; }
    .stDownloadButton > button { width: 100%; }
    .metric-card {
        background: #f0f4ff;
        border-radius: 8px;
        padding: 0.8rem 1rem;
        text-align: center;
    }
</style>
""",
    unsafe_allow_html=True,
)

# ── sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.image(
        "https://upload.wikimedia.org/wikipedia/commons/thumb/a/a7/"
        "Camponotus_flavomarginatus_ant.jpg/240px-Camponotus_flavomarginatus_ant.jpg",
        use_column_width=False,
        width=50,
    )
    st.title("🧬 Protein Domain Toolkit")
    st.caption("Dr. Arli Aditya Parikesit")
    st.caption("Department of Bioinformatics, i3L University")
    st.divider()

    tool = st.radio(
        "**Select Tool**",
        options=[
            "⚡ Full S. cerevisiae Workflow",
            "🔬 UniProt Mining",
            "📊 Domain Co-occurrence",
            "📈 XMGrace Output",
            "🧫 BED → FASTA",
            "🧹 Remove Redundant FASTA",
            "🔗 Network Duplicate Detection",
            "🗂️ Overlap Region Collapse",
        ],
    )
    st.divider()
    st.markdown("**References**")
    st.caption(
        "Parikesit et al. (2014). *Malaysian J. Fund. Appl. Sci.*, 10(2), 65–75."
        " [DOI](https://doi.org/10.11113/mjfas.v10n2.57)"
    )
    st.caption(
        "Parikesit et al. (2011). *Genes*, 2(4), 912–924."
        " [DOI](https://doi.org/10.3390/genes2040912)"
    )
    st.divider()
    st.caption("PhD dissertation (Leipzig 2012):")
    st.caption(
        "[figshare.964089](https://doi.org/10.6084/M9.FIGSHARE.964089)"
    )


# ════════════════════════════════════════════════════════════════════════════
#  Helper utilities
# ════════════════════════════════════════════════════════════════════════════

def _df_download_button(
    df: pd.DataFrame,
    label: str,
    filename: str,
    sep: str = "\t",
    mime: str = "text/plain",
) -> None:
    buf = io.BytesIO()
    df.to_csv(buf, sep=sep, index=False)
    st.download_button(label, data=buf.getvalue(), file_name=filename, mime=mime)


def _bytes_download_button(
    data: bytes | str,
    label: str,
    filename: str,
    mime: str = "text/plain",
) -> None:
    if isinstance(data, str):
        data = data.encode()
    st.download_button(label, data=data, file_name=filename, mime=mime)


def _load_df_from_upload(uploaded, sep_hint: str = "\t") -> pd.DataFrame | None:
    if uploaded is None:
        return None
    content = uploaded.read().decode("utf-8", errors="replace")
    sep = sep_hint
    if uploaded.name.endswith(".csv"):
        sep = ","
    try:
        return pd.read_csv(io.StringIO(content), sep=sep)
    except Exception as exc:
        st.error(f"Could not parse file: {exc}")
        return None


# ════════════════════════════════════════════════════════════════════════════
#  TOOL 0 — Full S. cerevisiae workflow  (recommended entry-point)
# ════════════════════════════════════════════════════════════════════════════

if tool == "⚡ Full S. cerevisiae Workflow":
    st.header("Full S. cerevisiae Proteomics Domain Workflow")
    st.markdown(
        """
        End-to-end demonstration using *Saccharomyces cerevisiae* (Baker's yeast,
        TaxID 559292) as a canonical model eukaryote.  Yeast is ideal because its
        proteome is exhaustively curated in UniProt Swiss-Prot, making it a
        dependable reference for domain distribution studies.

        **Pipeline steps**

        | Step | Description |
        |------|-------------|
        | 1 | Mine reviewed proteins + domain annotations from UniProt REST API |
        | 2 | Parse `ft_domain` fields; count pairwise domain co-occurrences |
        | 3 | Generate XMGrace-compatible `.agr` bar-chart file |
        | 4 | Render publication-quality heatmap (PNG) |
        """
    )

    col_left, col_right = st.columns(2)
    with col_left:
        max_proteins = st.slider("Max proteins to retrieve", 100, 1000, 300, 50,
                                 help="UniProt REST returns up to 500 per page.")
        top_n_heat   = st.slider("Top N domains in heatmap", 10, 40, 25)
    with col_right:
        top_n_agr = st.slider("Top N pairs in XMGrace chart", 5, 30, 20)
        agr_title = st.text_input(
            "XMGrace chart title",
            "Domain Co-occurrence in Saccharomyces cerevisiae",
        )

    if st.button("▶  Run Full Workflow", type="primary", use_container_width=True):
        progress = st.progress(0, text="Initialising …")
        log      = st.empty()

        try:
            # ── step 1 ──────────────────────────────────────────────────────
            log.info("Step 1/4 — Mining UniProt for S. cerevisiae …")
            progress.progress(10, text="Mining UniProt …")
            df = mine_uniprot(SCEREVISIAE_TAXID, reviewed=True,
                              max_results=max_proteins)
            st.session_state["uniprot_df"] = df
            progress.progress(30, text="UniProt data retrieved.")

            # stats row
            domain_col = next(
                (c for c in df.columns if "domain" in c.lower()), None
            )
            n_with_domains = 0
            if domain_col:
                n_with_domains = df[domain_col].notna().sum()
                n_with_domains -= (df[domain_col] == "").sum()

            m1, m2, m3 = st.columns(3)
            m1.metric("Proteins retrieved",          len(df))
            m2.metric("With domain annotations",     n_with_domains)
            m3.metric("Organism",                    "S. cerevisiae")

            with st.expander("Proteome preview (first 10 rows)"):
                st.dataframe(df.head(10), use_container_width=True)

            # ── step 2 ──────────────────────────────────────────────────────
            log.info("Step 2/4 — Computing domain co-occurrence …")
            progress.progress(50, text="Computing co-occurrence matrix …")
            matrix, cooc = compute_domain_cooccurrence(df)
            st.session_state["cooccurrence_matrix"]  = matrix
            st.session_state["cooccurrence_counter"] = cooc

            progress.progress(65)

            if matrix.empty:
                st.warning(
                    "No domain annotations found in the retrieved proteins. "
                    "Try increasing **Max proteins** or check the UniProt "
                    "connection."
                )
            else:
                st.subheader("Top Co-occurring Domain Pairs")
                top_pairs_data = [
                    {
                        "Rank":           i + 1,
                        "Domain A":       p[0],
                        "Domain B":       p[1],
                        "Co-occurrence":  c,
                    }
                    for i, (p, c) in enumerate(
                        [(pp, cc) for pp, cc in cooc.most_common(30)
                         if pp[0] != pp[1]][:25]
                    )
                ]
                if top_pairs_data:
                    st.dataframe(
                        pd.DataFrame(top_pairs_data),
                        use_container_width=True,
                        hide_index=True,
                    )

                # ── step 3 ──────────────────────────────────────────────────
                log.info("Step 3/4 — Generating XMGrace .agr file …")
                progress.progress(75, text="Generating XMGrace output …")
                agr_path = RESULTS_DIR / "scerevisiae_domain_cooccurrence.agr"
                agr_content = generate_xmgrace_cooccurrence(
                    cooc, agr_path, top_n=top_n_agr, title=agr_title
                )

                with st.expander("XMGrace .agr file preview (first 70 lines)"):
                    st.code(
                        "\n".join(agr_content.split("\n")[:70]),
                        language="bash",
                    )

                # ── step 4 ──────────────────────────────────────────────────
                log.info("Step 4/4 — Rendering heatmap …")
                progress.progress(90, text="Rendering heatmap …")
                hm_path = RESULTS_DIR / "scerevisiae_domain_heatmap.png"
                plot_cooccurrence_heatmap(
                    matrix, hm_path,
                    title=(
                        "Domain Co-occurrence — Saccharomyces cerevisiae "
                        "(UniProt Swiss-Prot reviewed)"
                    ),
                    top_n=top_n_heat,
                )
                progress.progress(100, text="Done.")
                log.success("Workflow complete.")

                st.subheader("Domain Co-occurrence Heatmap")
                st.image(str(hm_path), use_column_width=True)

                # ── downloads ────────────────────────────────────────────────
                st.subheader("Download Results")
                d1, d2, d3, d4 = st.columns(4)
                with d1:
                    _df_download_button(
                        df, "Proteome TSV",
                        "scerevisiae_proteome.tsv"
                    )
                with d2:
                    _df_download_button(
                        matrix, "Co-occurrence Matrix",
                        "cooccurrence_matrix.tsv"
                    )
                with d3:
                    with open(hm_path, "rb") as fh:
                        _bytes_download_button(
                            fh.read(), "Heatmap PNG",
                            "domain_cooccurrence_heatmap.png",
                            mime="image/png",
                        )
                with d4:
                    _bytes_download_button(
                        agr_content, "XMGrace .agr",
                        "scerevisiae_domain_cooccurrence.agr"
                    )

                st.info(
                    "**Render XMGrace file:**  \n"
                    "`xmgrace scerevisiae_domain_cooccurrence.agr`  \n"
                    "**Non-interactive PNG batch:**  \n"
                    "`gracebat scerevisiae_domain_cooccurrence.agr "
                    "-hdevice PNG -hardcopy -printfile output.png`"
                )

        except Exception as exc:
            log.error(f"Workflow failed: {exc}")
            st.error(f"Error: {exc}")
            st.exception(exc)


# ════════════════════════════════════════════════════════════════════════════
#  TOOL 1 — UniProt mining
# ════════════════════════════════════════════════════════════════════════════

elif tool == "🔬 UniProt Mining":
    st.header("UniProt Proteomics Data Mining")
    st.markdown(
        "Query the UniProt REST API for any organism and retrieve protein "
        "entries with domain annotations.  Default: *S. cerevisiae* (TaxID 559292)."
    )

    c1, c2 = st.columns(2)
    with c1:
        org_id = st.number_input(
            "Organism Taxonomy ID",
            min_value=1, value=SCEREVISIAE_TAXID,
            help="NCBI Taxonomy ID.  S. cerevisiae = 559292, H. sapiens = 9606"
        )
    with c2:
        max_r = st.slider("Max proteins", 50, 1000, 300, 50)

    reviewed = st.checkbox("Swiss-Prot reviewed only", value=True)

    if st.button("Mine UniProt", type="primary"):
        with st.spinner("Querying UniProt REST API …"):
            try:
                df = mine_uniprot(
                    organism_id=org_id,
                    reviewed=reviewed,
                    max_results=max_r,
                )
                st.session_state["uniprot_df"] = df
                st.success(f"Retrieved {len(df)} proteins.")

                st.dataframe(df.head(20), use_container_width=True)

                buf = io.BytesIO()
                df.to_csv(buf, sep="\t", index=False)
                st.download_button(
                    "Download TSV",
                    data=buf.getvalue(),
                    file_name=f"uniprot_{org_id}.tsv",
                    mime="text/plain",
                )
            except Exception as exc:
                st.error(f"Error: {exc}")
                st.exception(exc)


# ════════════════════════════════════════════════════════════════════════════
#  TOOL 2 — Domain co-occurrence
# ════════════════════════════════════════════════════════════════════════════

elif tool == "📊 Domain Co-occurrence":
    st.header("Domain Co-occurrence Analysis")
    st.markdown(
        "Compute a symmetric co-occurrence matrix for all protein domain pairs. "
        "Upload a UniProt TSV or use data from the mining step."
    )

    uploaded = st.file_uploader(
        "Upload UniProt TSV (or use previously mined data below)",
        type=["tsv", "txt", "csv"],
    )

    df: pd.DataFrame | None = None
    if uploaded:
        df = _load_df_from_upload(uploaded)
        if df is not None:
            st.info(f"Loaded {len(df)} rows, {len(df.columns)} columns.")
    elif "uniprot_df" in st.session_state:
        df = st.session_state["uniprot_df"]
        st.info(f"Using session data ({len(df)} proteins).")

    if df is not None:
        top_n = st.slider("Top N domains in heatmap", 5, 40, 20)

        if st.button("Compute Co-occurrence", type="primary"):
            with st.spinner("Computing …"):
                try:
                    matrix, cooc = compute_domain_cooccurrence(df)
                    st.session_state["cooccurrence_matrix"]  = matrix
                    st.session_state["cooccurrence_counter"] = cooc

                    if matrix.empty:
                        st.warning(
                            "No domain annotations found.  "
                            "Check that the TSV contains an 'ft_domain' column."
                        )
                    else:
                        st.success(
                            f"Co-occurrence matrix: "
                            f"{matrix.shape[0]} × {matrix.shape[1]} domains"
                        )

                        # top pairs table
                        st.subheader("Top 20 Co-occurring Pairs")
                        tbl = [
                            {"Domain A": p[0], "Domain B": p[1],
                             "Co-occurrence": c}
                            for p, c in cooc.most_common(30)
                            if p[0] != p[1]
                        ][:20]
                        if tbl:
                            st.dataframe(
                                pd.DataFrame(tbl),
                                use_container_width=True,
                                hide_index=True,
                            )

                        # heatmap
                        st.subheader("Co-occurrence Heatmap")
                        hm_path = RESULTS_DIR / "streamlit_heatmap.png"
                        plot_cooccurrence_heatmap(
                            matrix, hm_path, top_n=top_n,
                            title=f"Domain Co-occurrence (top {top_n} domains)",
                        )
                        st.image(str(hm_path), use_column_width=True)

                        # downloads
                        c_a, c_b = st.columns(2)
                        with c_a:
                            mat_buf = io.BytesIO()
                            matrix.to_csv(mat_buf, sep="\t")
                            st.download_button(
                                "Download Matrix TSV",
                                data=mat_buf.getvalue(),
                                file_name="cooccurrence_matrix.tsv",
                                mime="text/plain",
                            )
                        with c_b:
                            with open(hm_path, "rb") as fh:
                                st.download_button(
                                    "Download Heatmap PNG",
                                    data=fh.read(),
                                    file_name="domain_cooccurrence_heatmap.png",
                                    mime="image/png",
                                )

                except Exception as exc:
                    st.error(f"Error: {exc}")
                    st.exception(exc)
    else:
        st.info("Mine data first (UniProt Mining tab) or upload a TSV file.")


# ════════════════════════════════════════════════════════════════════════════
#  TOOL 3 — XMGrace output
# ════════════════════════════════════════════════════════════════════════════

elif tool == "📈 XMGrace Output":
    st.header("XMGrace .agr Bar-Chart Generator")
    st.markdown(
        "Produce an XMGrace project file (`.agr`) showing the top-N co-occurring "
        "domain pairs as a colour-filled bar chart.  XMGrace is widely used in "
        "structural biology and computational genomics for publication-quality plots."
    )

    # try session state first
    cooc: Counter | None = st.session_state.get("cooccurrence_counter")
    matrix: pd.DataFrame | None = st.session_state.get("cooccurrence_matrix")

    # allow manual upload
    uploaded_mat = st.file_uploader(
        "Or upload co-occurrence matrix TSV (optional)",
        type=["tsv", "txt"],
    )
    if uploaded_mat:
        matrix = pd.read_csv(
            io.StringIO(uploaded_mat.read().decode()), sep="\t", index_col=0
        )
        cooc = Counter()
        for d1 in matrix.index:
            for d2 in matrix.columns:
                cnt = int(matrix.loc[d1, d2])
                if cnt > 0 and str(d1) <= str(d2):
                    cooc[(str(d1), str(d2))] = cnt

    if cooc is None and matrix is not None:
        cooc = Counter()
        for d1 in matrix.index:
            for d2 in matrix.columns:
                cnt = int(matrix.loc[d1, d2])
                if cnt > 0 and str(d1) <= str(d2):
                    cooc[(str(d1), str(d2))] = cnt

    c1, c2 = st.columns(2)
    with c1:
        top_n = st.slider("Top N domain pairs", 5, 30, 20)
    with c2:
        plot_title = st.text_input(
            "Plot title",
            "Domain Co-occurrence in Saccharomyces cerevisiae",
        )

    if cooc:
        if st.button("Generate XMGrace .agr", type="primary"):
            agr_path = RESULTS_DIR / "domain_cooccurrence.agr"
            content  = generate_xmgrace_cooccurrence(
                cooc, agr_path, top_n=top_n, title=plot_title
            )
            st.success("XMGrace file generated.")

            with st.expander("File preview (first 80 lines)"):
                st.code(
                    "\n".join(content.split("\n")[:80]), language="bash"
                )

            _bytes_download_button(
                content,
                "Download .agr File",
                "domain_cooccurrence.agr",
            )

            st.info(
                "**Interactive:**  `xmgrace domain_cooccurrence.agr`\n\n"
                "**Non-interactive PNG export:**\n"
                "```\ngracebat domain_cooccurrence.agr \\\n"
                "  -hdevice PNG -hardcopy -printfile output.png\n```"
            )
    else:
        st.info(
            "Run the **Domain Co-occurrence** step first, or upload a "
            "co-occurrence matrix TSV above."
        )


# ════════════════════════════════════════════════════════════════════════════
#  TOOL 4 — BED → FASTA
# ════════════════════════════════════════════════════════════════════════════

elif tool == "🧫 BED → FASTA":
    st.header("BED Coordinates → FASTA Extraction")
    st.markdown(
        "Convert a BED file of genomic coordinates to FASTA sequences. "
        "Requires a local genome directory containing per-chromosome FASTA "
        "files (`<chrom>.fa`).  Uses **pyfaidx** when available; otherwise "
        "falls back to the legacy `fastacmd` binary."
    )

    bed_file  = st.file_uploader("Upload BED file", type=["bed", "txt"])
    genome_dir = st.text_input(
        "Genome directory (server-side path)",
        placeholder="/path/to/genome/chromosomes/",
    )

    if bed_file and genome_dir:
        if st.button("Convert BED → FASTA", type="primary"):
            with tempfile.NamedTemporaryFile(
                suffix=".bed", delete=False, mode="w"
            ) as tmp:
                tmp.write(bed_file.getvalue().decode())
                tmp_path = tmp.name

            out_path = RESULTS_DIR / "bed_output.fasta"
            try:
                n = bed_to_fasta(tmp_path, genome_dir, str(out_path))
                st.success(f"Written {n} FASTA records.")
                if out_path.exists():
                    with open(out_path, "rb") as fh:
                        _bytes_download_button(
                            fh.read(), "Download FASTA",
                            "bed_to_fasta.fa", mime="text/plain"
                        )
            except Exception as exc:
                st.error(f"Error: {exc}")
            finally:
                os.unlink(tmp_path)
    elif not genome_dir:
        st.info("Provide the server-side genome directory path.")


# ════════════════════════════════════════════════════════════════════════════
#  TOOL 5 — Remove redundant FASTA
# ════════════════════════════════════════════════════════════════════════════

elif tool == "🧹 Remove Redundant FASTA":
    st.header("Remove Redundant FASTA Sequences")
    st.markdown(
        "Collapse identical sequences.  When multiple entries share the same "
        "residue string, their identifiers are concatenated with `;` and a "
        "single representative record is emitted — reproducing the behaviour "
        "of the original Perl script."
    )

    fasta_file = st.file_uploader(
        "Upload FASTA file",
        type=["fa", "fasta", "fna", "faa", "txt"],
    )

    if fasta_file:
        if st.button("Remove Redundancy", type="primary"):
            with tempfile.NamedTemporaryFile(
                suffix=".fa", delete=False, mode="wb"
            ) as tmp_in:
                tmp_in.write(fasta_file.read())
                tmp_in_path = tmp_in.name

            tmp_out_path = tmp_in_path + "_nr.fa"
            try:
                total, unique = remove_redundant_fasta(
                    tmp_in_path, tmp_out_path
                )
                c1, c2, c3 = st.columns(3)
                c1.metric("Input sequences",   total)
                c2.metric("Unique sequences",  unique)
                c3.metric("Removed",           total - unique)

                with open(tmp_out_path, "rb") as fh:
                    _bytes_download_button(
                        fh.read(),
                        "Download Non-redundant FASTA",
                        "nonredundant.fa",
                    )
            except Exception as exc:
                st.error(f"Error: {exc}")
            finally:
                for p in [tmp_in_path, tmp_out_path]:
                    if os.path.exists(p):
                        os.unlink(p)


# ════════════════════════════════════════════════════════════════════════════
#  TOOL 6 — Network duplicate detection
# ════════════════════════════════════════════════════════════════════════════

elif tool == "🔗 Network Duplicate Detection":
    st.header("Duplicate Edge Detection in Biological Networks")
    st.markdown(
        "Identify duplicate edges in a tab-separated network file. "
        "Edge direction is normalised (A→B ≡ B→A) before counting, "
        "reproducing the awk pipeline from `Duplicates in Biological Network.sh`."
    )

    net_file = st.file_uploader(
        "Upload network file (TSV)", type=["tsv", "txt", "tab"]
    )

    c1, c2, c3 = st.columns(3)
    with c1:
        col_a = st.number_input("Node A column index",       0, 20, 0)
    with c2:
        col_b = st.number_input("Interaction column index",  0, 20, 1)
    with c3:
        col_c = st.number_input("Node B column index",       0, 20, 5)

    if net_file:
        if st.button("Detect Duplicates", type="primary"):
            with tempfile.NamedTemporaryFile(
                suffix=".txt", delete=False, mode="wb"
            ) as tmp:
                tmp.write(net_file.read())
                tmp_path = tmp.name

            out_path = RESULTS_DIR / "deduped_network.tsv"
            try:
                result_df = detect_network_duplicates(
                    tmp_path, str(out_path),
                    col1=int(col_a), col2=int(col_b), col3=int(col_c),
                )
                st.success(f"Processed — {len(result_df)} unique edges.")
                st.dataframe(result_df.head(20), use_container_width=True)
                _df_download_button(
                    result_df, "Download Deduplicated Network",
                    "deduped_network.tsv"
                )
            except Exception as exc:
                st.error(f"Error: {exc}")
            finally:
                os.unlink(tmp_path)


# ════════════════════════════════════════════════════════════════════════════
#  TOOL 7 — Overlap region collapse
# ════════════════════════════════════════════════════════════════════════════

elif tool == "🗂️ Overlap Region Collapse":
    st.header("Overlapping Domain Region Collapse")
    st.markdown(
        "Merge overlapping BLAST/domain-hit intervals per protein and retain "
        "the best hit (lowest e-value) from each cluster — the core algorithm "
        "from `Domain cooccurence.sh` and `Domain cooccurence-o.sh`."
    )

    overlap_file = st.file_uploader(
        "Upload domain overlap list (tab-separated)",
        type=["tsv", "txt", "list"],
    )

    if overlap_file:
        if st.button("Collapse Overlapping Regions", type="primary"):
            with tempfile.NamedTemporaryFile(
                suffix=".txt", delete=False, mode="wb"
            ) as tmp:
                tmp.write(overlap_file.read())
                tmp_path = tmp.name

            try:
                result_df = process_overlap_regions(
                    tmp_path,
                    str(RESULTS_DIR / "min_nmbr.txt"),
                    str(RESULTS_DIR / "max_nmbr.txt"),
                    str(RESULTS_DIR / "region_limit.txt"),
                )
                st.success(
                    f"Collapsed to {len(result_df)} non-overlapping regions."
                )
                st.dataframe(result_df.head(20), use_container_width=True)

                c1, c2 = st.columns(2)
                with c1:
                    _df_download_button(
                        result_df, "Download Processed Regions",
                        "processed_regions.tsv"
                    )
                with c2:
                    if (RESULTS_DIR / "region_limit.txt").exists():
                        with open(RESULTS_DIR / "region_limit.txt", "rb") as fh:
                            _bytes_download_button(
                                fh.read(), "Download Region Limits",
                                "region_limit.tsv"
                            )
            except Exception as exc:
                st.error(f"Error: {exc}")
            finally:
                os.unlink(tmp_path)


# ── footer ────────────────────────────────────────────────────────────────────
st.divider()
st.caption(
    "Protein Domain Analysis Toolkit v2.0  |  "
    "Dr. Arli Aditya Parikesit  |  i3L University Jakarta  |  "
    "ORCID: [0000-0001-8716-3926](https://orcid.org/0000-0001-8716-3926)  |  "
    "Dissertation: [figshare.964089](https://doi.org/10.6084/M9.FIGSHARE.964089)"
)
