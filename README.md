# Elakatotrichophyceae phylogenomics

Scripts and documented analysis steps for the MBE major revision concerning
Elakatotrichophyceae. The repository URL retains its historical spelling,
`Elakatophyceae-Phylogenomics`; the formal class name is **Elakatotrichophyceae**.

The revised phylogenetic dataset contains **68 taxa, 388 loci and 130,892
aligned amino-acid sites**. See [the revised workflow](WORKFLOW_REVISED_388LOCI.md)
for the actual selection, alignment, concatenation and tree-inference stages.

## Data and code

Final alignments, supermatrix, taxon/locus lists, partitions and species trees
are deposited at Zenodo: https://doi.org/10.5281/zenodo.23855065.
The record is published with restricted files; reviewer access is provided
through a confidential sharing link.

The deposit contains two data ZIPs, a README, an upload list and a checksum
file. `TREE_INDEX.xlsx` maps the four main species trees to Figure 1/S-PHY2 and
outgroup-comparison trees to Figure 2A.

## Script versions and execution order

The root-level scripts are retained original-submission workflow material.
They are not an automated implementation of the revised 388-locus analysis.
Versions are not pinned in those scripts; executable names are not version records.

| Order / script | Actual input and operation | Main output / dependencies |
| --- | --- | --- |
| 1 `step1.sh` | Existing OrthoFinder `Orthogroups.tsv`; select 5-of-5 and 4-of-5 single-copy groups | `SOG_5of5_1copy.txt`, `SOG_4of5_1copy.txt`; Python, pandas |
| 2 `step2.sh` | Listed `Orthogroup_Sequences/*.fa`; align seeds and build/press HMM profiles | `SOG468_aln/`, `SOG468_hmm/`, `SOG468.hmmdb`; Bash, MAFFT, HMMER |
| 3 `step3.sh` | HMM library and sample `.pep` files; `hmmsearch -E 1e-6 --domE 1e-6` | `hmmsearch_tbl/`, `hmmsearch_domtbl/`, logs; HMMER |
| 4 `step4.sh` | Domain tables and proteins; one best hit per sample/locus by i-Evalue, then domain score | `best_hits.tsv`, `step4_loci_fa/`, `per_locus_counts.tsv`; Python |
| 5 `step5.sh` | Per-locus FASTAs; retain at least 54 sequences | `step5_loci_ge80/`; Bash, grep |
| 6 `step6.sh` | Selected FASTAs; MAFFT `--auto`, trimAl `-gt 0.70 -cons 50` | `step6_aln/`, `step6_trim/`; MAFFT, trimAl |
| Between 6 and 7 `rename.sh` | Trimmed headers `sample|sequence`; retain the sample part | `step6_trim_renamed/`; awk |
| 7 `step7.sh` | Renamed alignments; IQ-TREE ModelFinder and 1,000 UFBoot replicates, TreeShrink q=0.05, occupancy/length filter | Gene trees, TreeShrink outputs, final alignments and filter TSV; IQ-TREE 2, TreeShrink, Python |
| Separate `eleka.ASR.R` | `characters.txt`, `ASR.tre`; `chronos(lambda=0.1)`, discrete `ace` ER/SYM/ARD selection by AIC | `ASR.pdf` and in-memory model results; R, ape, phytools, plotrix |

These are server-specific script fragments, not a command to run sequentially
without preparation. Existing inputs and directories are required. Paths differ
between scripts (including `OrthoFinder` versus `Orthofinder`); `step7.sh` has
`WORKDIR="cd /mnt/..."`, which must be corrected in a separate working copy
before use. That script also clears output directories. `step4.sh` appends FASTA
records, so rerunning it in a populated output directory can duplicate records.
`eleka.ASR.R` implements the original-submission ASR workflow.

## Revised analysis scripts

`revised_scripts/` contains nine scripts for the revised workflow.
`SOURCE_MANIFEST.tsv` lists their source paths and SHA-256 values. The revised guide identifies which commands and
parameters are supported by run records, and which upstream utilities have
only their source code available. In particular, the five-seed tree-guided
selector supersedes the root-level `step4.sh` single-best-hit selection.

Linux/Bash, Python with Biopython, MAFFT, trimAl, HMMER, IQ-TREE and Java/ASTRAL
are needed for the corresponding stages. IQ-TREE **2.4.0** is verified in the
final C60 run log; ASTRAL **5.7.8** is the main unconstrained analysis and
**5.6.9** the separate constrained analysis. Exact historical versions of the
other programs are not established by the inspected scripts, and are not
assigned here. OrthoFinder/DIAMOND outputs precede this workflow; the original
protein-prediction and orthogroup-generation commands are not supplied here.

This repository documents the recorded workflow and provides reusable scripts.
The minimal deposit supplies final phylogenetic inputs and results, rather
than all upstream HMM databases and candidate-search intermediates. No claim
is made that the legacy scripts provide a one-command reconstruction of all
revised analyses. The workflow guide covers the phylogenomic analyses.
