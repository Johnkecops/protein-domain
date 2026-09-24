# Domain Co-occurrence — Top 20 Pairs (S. cerevisiae, UniProt Swiss-Prot reviewed)

Ground truth, extracted verbatim from the `# pair N: ...` comment block that
`generate_xmgrace_cooccurrence()` writes at the end of
[`results/domain_cooccurrence.agr`](domain_cooccurrence.agr). `n` is the exact
protein count from `compute_domain_cooccurrence()` — not a pixel-color read of
a heatmap PNG.

Use this table to answer "does domain A co-occur with domain B" questions.
The heatmap images (`scerevisiae_domain_heatmap.png`, `streamlit_heatmap.png`)
are for visual browsing only: an AI or human reading color intensity off those
PNGs cannot reliably distinguish n=1 from n=0, which produced several
unsupported "domain X relates to domain Y" claims in `graphify-out/graph.json`.
Any such claim not listed below is a heatmap-reading artifact, not a real
co-occurrence.

| Rank | Domain A | Domain B | n |
|---|---|---|---|
| 1 | Helicase ATP-binding | Helicase C-terminal | 7 |
| 2 | AGC-kinase C-terminal | Protein kinase | 3 |
| 3 | Myosin N-terminal SH3-like | Myosin motor | 2 |
| 4 | Integrase catalytic | RNase H Ty1/copia-type | 2 |
| 5 | Integrase catalytic | Reverse transcriptase Ty1/copia-type | 2 |
| 6 | RNase H Ty1/copia-type | Reverse transcriptase Ty1/copia-type | 2 |
| 7 | ACT-like 1 | ACT-like 2 | 1 |
| 8 | ATP-grasp 1 | ATP-grasp 2 | 1 |
| 9 | ATP-grasp 1 | MGS-like | 1 |
| 10 | ATP-grasp 2 | MGS-like | 1 |
| 11 | Protein kinase | UBA | 1 |
| 12 | N-terminal Ras-GEF | Ras-GEF | 1 |
| 13 | IQ | Myosin N-terminal SH3-like | 1 |
| 14 | IQ | Myosin motor | 1 |
| 15 | Inhibitor I9 | Peptidase S8 | 1 |
| 16 | ATP-grasp | Biotin carboxylation | 1 |
| 17 | ATP-grasp | Biotinyl-binding | 1 |
| 18 | ATP-grasp | Pyruvate carboxyltransferase | 1 |
| 19 | Biotin carboxylation | Biotinyl-binding | 1 |
| 20 | Biotin carboxylation | Pyruvate carboxyltransferase | 1 |

**Regenerate after re-running the pipeline:** re-run
`_run_scerevisiae_workflow()` (or `python protein_domain_toolkit.py`), then
copy the new `# pair N: ...` block from `domain_cooccurrence.agr` into the
table above — the workflow does not write this file automatically.

## Resolving the graph's AMBIGUOUS domain-pair questions

Cross-checked against the table above:

- **Helicase C-terminal ↔ SANT 2** — SANT 2 does not appear in the top-20 at
  all. No supported relationship; the graph edge is a heatmap-reading
  artifact and should be disregarded or removed.
- **Integrase catalytic ↔ Helicase ATP-binding domain 1/2** — Integrase
  catalytic co-occurs only with RNase H Ty1/copia-type and Reverse
  transcriptase Ty1/copia-type (both n=2, rows 4–5), never with Helicase
  ATP-binding. The two domain families sit in unrelated top-20 rows; treat
  the graph edge as unsupported.
- **EF-hand ↔ Pyruvate carboxyltransferase** — EF-hand does not appear in the
  top-20 at all. Pyruvate carboxyltransferase's only real partners are
  ATP-grasp and Biotin carboxylation (row 18, 20). Unsupported.

All three unsupported edges share a cause: the vision extraction read
low-intensity heatmap cells as plausible-looking co-occurrences without
access to the underlying counts. The counts were available the whole time,
in this repo, at the bottom of `domain_cooccurrence.agr`.
