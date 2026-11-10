# GRIFO

**GRIFO** is a Jupyter-notebook pipeline for **Oxford Nanopore amplicon sequencing** (eDNA
metabarcoding). It takes demultiplexed reads, filters them, clusters them into consensus
sequences without a reference, pools those into OTUs shared across samples, assigns their
taxonomy with BLAST, and draws an interactive haplotype network for every OTU.

```
raw reads ──► 1 QC (NanoFilt) ──► 3/4 per-sample clustering (minimap2 + spoa)
          ──► 5 OTUs across samples (vsearch) + taxonomy (BLAST, haplonet.py)
          ──► 6 haplotype networks (MAFFT, haplonet.py) ──► one HTML page
```

---

## Contents

| file | what it is |
|---|---|
| `grifo.ipynb` | the pipeline, one section per step; run it top to bottom |
| `config.ini` | paths and settings read by the notebook |
| `haplonet.py` | BLAST taxonomy and haplotype networks; used by the notebook, also runs on its own |
| `clustering_tools.py` | clustering code; **written by the notebook** (section 4), don't edit it by hand |
| `cluster.py` | the same clustering code as a standalone command-line script |

Folders you provide or the pipeline creates (keep them out of git, see [What to commit](#what-to-commit)):

| folder / file | contents |
|---|---|
| `data/` | your demultiplexed reads (input) |
| `MIDORI2_…fasta` | reference database (input) |
| `results/qc/` | filtered reads, clusters and read counts per sample, MultiQC report |
| `results/taxonomy/` | OTU table, OTU sequences, BLAST hits, BLAST database |
| `results/haplotype_networks/` | `haplotype_networks.html` and the alignments per OTU |

---

## Installation

GRIFO needs Python with numpy, pandas and scikit-learn, Jupyter, and these command-line tools on
your `PATH`: NanoFilt, NanoStat, MultiQC, minimap2, spoa, vsearch, BLAST+ and MAFFT. The easiest
way is one conda environment (Miniforge or Mambaforge):

```bash
mamba create -n grifo -c conda-forge -c bioconda \
    python=3.11 numpy pandas scikit-learn jupyterlab \
    nanofilt nanostat multiqc minimap2 spoa vsearch blast mafft
mamba activate grifo
```

Check that everything is found:

```bash
for t in NanoFilt NanoStat multiqc minimap2 spoa vsearch blastn makeblastdb mafft; do
    printf '%-12s ' $t; command -v $t || echo MISSING
done
```

The haplotype-network page loads its drawing library (vis-network) from the internet, so you need
a connection to *view* it, not to make it.

### Reference database

Taxonomy is assigned against a reference FASTA. GRIFO was built around
[MIDORI2](https://www.reference-midori.info/) (mitochondrial references from GenBank), e.g. the
16S file `MIDORI2_TOTAL_NUC_GB273_lrRNA_RAW.fasta`. Download the marker that matches your primers,
unzip it, and put it next to the notebook.

`haplonet.py` also reads SINTAX/UTAX, QIIME/UNITE, SILVA/PR2 and plain `>ID Genus species` headers.

---

## Tutorial

This walks through one complete run. Expect the first run to take longest: the BLAST database is
built once from the reference FASTA and reused afterwards.

### 1. Put your data in place

Use one folder per sample under `data/`, with the gzipped FASTQ files directly inside. The folder
name becomes the sample name everywhere downstream:

```
v_0.7/
├── grifo.ipynb
├── haplonet.py
├── config.ini
├── MIDORI2_TOTAL_NUC_GB273_lrRNA_RAW.fasta
└── data/
    ├── barcode01/
    │   ├── fastq_runid_..._0_0.fastq.gz
    │   └── fastq_runid_..._1_0.fastq.gz
    └── barcode02/
        └── ...
```

This is the layout MinKNOW / Dorado produce after demultiplexing (`fastq_pass/barcodeXX/`).

> **Watch out:** every folder under `data/` that contains `*.fastq.gz` files becomes a sample,
> including `data/` itself if you leave FASTQ files at its top level, and every FASTQ file in a
> sample folder is used. Don't keep merged copies (e.g. `sample4.fastq.gz`) next to the original
> chunks, or those reads are counted twice.

### 2. Write `config.ini`

A minimal file:

```ini
[paths]
path_data = ./data
path_to_blastdb = MIDORI2_TOTAL_NUC_GB273_lrRNA_RAW.fasta

[NanoFilt]
minlength = 100      ; shortest read kept (bp)
maxlength = 800      ; longest read kept (bp)
qscore = 20          ; minimum mean read quality (whole number)

[clustering]
niter = 5            ; clustering iterations
```

Set `minlength` / `maxlength` around your expected amplicon length. All other settings have
defaults; the full list is in [Settings](#settings).

### 3. Run the notebook

```bash
mamba activate grifo
jupyter lab grifo.ipynb
```

Run the cells from top to bottom (**Kernel → Restart & Run All**), or section by section as below.

#### Section 1 – Setup
Reads `config.ini`, prints its sections and checks that `data/` exists. Check that the printed
settings are what you expect.

#### Section 2 – Read QC
Filters every FASTQ chunk with NanoFilt, merges each sample's chunks, and writes a NanoStat +
MultiQC report. The cell ends with a table:

```
sample                     raw  filtered    kept
barcode01                17999      9475   52.6%
barcode02               435086    313833   72.1%
```

Open `results/qc/multiqc/multiqc_report.html` to look at read lengths and quality. Samples that
keep far fewer reads than the others usually have a length or quality problem; adjust
`[NanoFilt]` and rerun this section.

#### Section 3 – Clustering input
Converts each sample's merged reads to `results/qc/<sample>/<sample>_concatenated.csv`.
From here on, everything works from the files in `results/qc`, so after a restart you can skip
section 2.

#### Section 4 – Clustering
Clusters each sample's reads into consensus sequences (reference-free: minimap2 distances,
OPTICS clustering, spoa consensus), then assigns every read to its closest cluster. Several
samples run in parallel (`PARALLEL_SAMPLES` and `THREADS_PER_SAMPLE` at the top of the cell). This is
the slowest step. The result:

```
sample                clusters     reads  assigned
barcode01                   29      9475      9427  (99.5%)
```

`assigned` should be close to 100%. Per sample, `results/qc/<sample>/` now holds
`_clusters.csv`, `_read_assignments.csv` and `_cluster_counts.csv`.

#### Section 5 – Taxonomy and OTU table
1. The clusters of all samples are pooled and clustered again at 97% with vsearch, so the same
   organism gets the same OTU in every sample. OTUs are numbered by total reads (`OTU_1` is the
   most abundant).
2. Each OTU is BLASTed against the reference (`haplonet.py`). If the best hit leaves ≥ 10 bases at
   an end unmatched (a primer concatemer or chimeric piece), the OTU is trimmed to the matched part.
3. The taxonomy is the deepest rank on which ≥ 60% of the best hits agree. Ranks the identity
   doesn't support are then cleared (species only at ≥ 98%, genus ≥ 95%, family ≥ 90%,
   order ≥ 85%).

The second code cell prints three parsed reference headers. **Check them before the search
runs**: every rank should show a sensible name. Empty ranks mean the header format isn't
recognised.

The last cell prints OTUs per assigned rank and shows the first rows of the table. The main result
is `results/taxonomy/otu_table.csv`:

| column | meaning |
|---|---|
| `barcode01` … | reads per sample |
| `total_reads` | reads over all samples |
| `domain` … `species` | taxonomy after the identity cutoffs |
| `assigned_rank` | deepest rank assigned, or `unassigned` |
| `consensus` | the consensus before the identity cutoffs, e.g. `Pomatoschistus microps (species)` |
| `identity`, `query_cov` | identity (%) and query cover of the best hit |
| `n_top_hits`, `support` | number of hits at the best identity, and the fraction that agree |
| `taxa_best_hits` | taxa of those best hits with their counts |
| `best_hit` | accession of the best hit |
| `query_region`, `trimmed_bp` | part of the OTU that was analysed, and how many bases were cut |
| `sequence` | the OTU sequence (untrimmed) |

**Reading the taxonomy.** A row with `genus = Pomatoschistus`, empty `species`, `identity = 97.3`
and `consensus = Pomatoschistus microps (species)` means: all best hits are *P. microps*, but at
97.3% identity the species is not accepted by the 98% cutoff. Whether to trust it depends on your
marker and how complete the reference is for your region. To keep species names below 98%, set
`rank_cutoffs = no` or lower `species_identity`.

#### Section 6 – Haplotype networks
For every OTU, its BLAST hits are aligned with the OTU (MAFFT), identical sequences are collapsed
into haplotypes, and a minimum spanning network is drawn. It reuses the hits from section 5, so
BLAST doesn't run again. Open

```
results/haplotype_networks/haplotype_networks.html
```

in a browser.

- **Pick an OTU** from the dropdown, or step through them with ← → (or ↑ ↓).
- **Filter** the OTU list by ID or taxon name, e.g. `Gobiidae`.
- The **red node** is your OTU. The other nodes are reference haplotypes, coloured by taxon;
  node size reflects the number of sequences, and the numbers on the edges are mutational steps.
- **Click a node** to list its accessions, identity to the OTU and full lineages.
- **Find** highlights the nodes matching a taxon or accession.
- **Relax / Freeze layout** and dragging nodes tidy up the picture; **Save PNG** exports it.

The side panel shows the OTU's read count, the OTU-table taxonomy and the consensus of the best
hits, plus a **Note** when the OTU was trimmed, had hits on N-runs discarded, or had no hits.

**What the network tells you.** An OTU that sits on (0 steps) or right next to the haplotypes of
one species is a confident identification. An OTU between two species, or far from all of them,
points to a species missing from the reference, a sequencing artefact, or closely related
species the marker can't separate.

### 4. Check a result at NCBI

To compare an OTU with NCBI's web BLAST, copy its `sequence` from `otu_table.csv` (or `otus.fasta`)
of the **same run**. OTU numbers are assigned by read count in each run, so `OTU_1` today can be a
different sequence from `OTU_1` in an earlier run or results folder.

Keep in mind that NCBI searches all of GenBank, while GRIFO only searches your reference file. An OTU
with *no hits* here may be off-target DNA (bacteria, nuclear DNA) that a mitochondrial database
cannot identify.

---

## Settings

All in `config.ini`; everything except `[paths]`, `[NanoFilt]` and `[clustering]` is optional.

```ini
[vsearch]
otu_identity = 0.97          ; identity at which sample clusters are pooled into OTUs
threads = 8

[haplonet]
db = MIDORI2_TOTAL_NUC_GB273_lrRNA_RAW.fasta   ; default: [paths] path_to_blastdb
blastdb =                    ; prebuilt BLAST database prefix; empty = build in results/taxonomy/blastdb
blast_task = megablast       ; dc-megablast or blastn find more divergent hits (slower)
threads = 8
top = 100                    ; hits kept per OTU
oversample = 5               ; extra hits requested, as some are discarded
min_id = 0.70                ; minimum hit identity
min_qcov = 80                ; minimum query cover of a hit (%)
evalue = 1e-5
max_ambig = 0.02             ; discard hits on reference regions with more N/IUPAC codes
trim_query = yes             ; trim OTUs to the part covered by their best hit
trim_min = 10                ; ...if at least this many bases are uncovered
tie_delta = 0.0              ; pool hits within this many %-points of the best
lca_cutoff = 0.6             ; fraction of best hits that must agree (1.0 = strict LCA)
rank_cutoffs = yes           ; clear ranks the identity doesn't support
species_identity = 0.98
genus_identity = 0.95
family_identity = 0.90
order_identity = 0.85
network_min_reads = 0        ; skip OTUs with fewer reads in the network page
mafft_mode = auto            ; auto, localpair, globalpair, genafpair
gaps = fifth                 ; gaps as fifth state, missing, or remove gap columns
max_label_taxa = 2           ; taxa named on a node before "+N more"
title = GRIFO haplotype networks
```

The BLAST database is rebuilt only if the reference FASTA is newer than it. If you built one
yourself, point `blastdb` to it (the name without `.nsq`/`.nin`); build it **without**
`-parse_seqids`, because MIDORI IDs aren't valid BLAST sequence IDs.

---

## Running `haplonet.py` on its own

`haplonet.py` works on any OTU table (CSV/TSV with an ID and a sequence column), a FASTA file or a
single sequence, without the notebook:

```bash
python haplonet.py -i results/taxonomy/otu_table.csv \
                   -d MIDORI2_TOTAL_NUC_GB273_lrRNA_RAW.fasta \
                   -o results_haplonet \
                   --blastdb results/taxonomy/blastdb/reference

python haplonet.py -i otus.fasta -d reference.fasta -o out --otus OTU_3,OTU_12
python haplonet.py --help        # all options
```

It writes the same HTML page, a `taxonomic_assignments.tsv` and a copy of the input table with
`hn_*` columns added.

---

## Troubleshooting

| problem | likely cause / fix |
|---|---|
| `Required tool(s) not found on PATH` | the conda environment isn't active in the Jupyter kernel; start Jupyter from the activated environment |
| a sample called `data` appears | FASTQ files at the top level of `data/`; move them into sample folders |
| reads counted twice | merged FASTQ files next to the original chunks in a sample folder |
| empty ranks in the header check (section 5) | reference header format not recognised; use a supported format or a `--taxonomy` TSV with `haplonet.py` |
| `BLAST hit IDs … are not in …` | `blastdb` was built from a different FASTA than `db` |
| many `unassigned` OTUs | wrong marker database, off-target amplification, or too strict cutoffs; look at `identity` and `consensus`, try `blast_task = dc-megablast` |
| species missing although `consensus` names one | the best identity is below `species_identity` |
| network page is blank | no internet connection (vis-network loads from a CDN) |
| results differ slightly between runs | clustering samples reads at random; numbers of clusters can vary a little |

---

## What to commit

Reads, references and results are large; keep them out of the repository. A suggested
`.gitignore`:

```gitignore
data/
results*/
*.zip
MIDORI2_*
*.udb
*.n??
__pycache__/
.ipynb_checkpoints/
.DS_Store
.~lock.*
```

Consider clearing the notebook outputs before committing (**Edit → Clear Outputs of All Cells**),
which keeps the diff small and avoids committing local paths.

---

## License

The clustering code (`cluster.py`, and `clustering_tools.py` written by the notebook) is a port of
code that is Copyright © 2020 Bilgenur Baloglu, Centre for Biodiversity Genomics, University of
Guelph, distributed under the
[GNU Lesser General Public License v3](https://www.gnu.org/licenses/lgpl-3.0.html), and is
distributed under the same license.

## Citing

If you use GRIFO, please also cite the tools it builds on: NanoFilt and NanoStat (NanoPack),
MultiQC, minimap2, spoa, vsearch, BLAST+, MAFFT, scikit-learn, and the reference database you used
(e.g. MIDORI2). The minimum spanning network follows Bandelt et al. (1999).
