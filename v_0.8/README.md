# GRIFO

**GRIFO** is a Jupyter-notebook pipeline for **Oxford Nanopore amplicon sequencing** (eDNA
metabarcoding). It takes demultiplexed reads, filters them, clusters them into consensus
sequences without a reference, pools those into OTUs shared across samples, assigns their
taxonomy with BLAST, draws an interactive haplotype network for every OTU, and packs the results
with your sample metadata into a Darwin Core Archive for GBIF and OBIS.

```
raw reads ──► 1 QC (NanoFilt) ──► 3/4 per-sample clustering (minimap2 + spoa)
          ──► 5 OTUs across samples (vsearch) + taxonomy (BLAST, haplonet.py)
          ──► 6 haplotype networks (MAFFT, haplonet.py) ──► one HTML page
          ──► 7 Darwin Core Archive (edna2dwca.py + sample metadata) ──► zip for GBIF/OBIS
          ──► 8 Krona plot of the taxonomy (ktImportText) ──► one HTML page
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
| `edna2dwca.py` | converter to a Darwin Core Archive, copied from [MarGenoDaT](https://github.com/pascalhabluetzel/MarGenoDaT) (section 7) |
| `metadata/sample_metadata_mock.csv` | **mock** sample metadata for testing section 7; replace it with your own |
| `metadata/manual_annotations.csv` | your manual taxonomic annotations, exported from the network page (section 6) and used in the archive (section 7); there only after your first export |
| `environment.yml` | the `grifo` conda environment |

Folders you provide or the pipeline creates (keep them out of git, see [What to commit](#what-to-commit)):

| folder / file | contents |
|---|---|
| `data/` | your demultiplexed reads (input) |
| `MIDORI2_…fasta` | reference database (input) |
| `results/qc/` | filtered reads, clusters and read counts per sample, MultiQC report |
| `results/taxonomy/` | OTU table, OTU sequences, BLAST hits, BLAST database |
| `results/haplotype_networks/` | `haplotype_networks.html` and the alignments per OTU |
| `results/dwca/` | Darwin Core Archive (`.zip`), its QC report and the converter's input files |
| `results/krona/` | `krona.html` and its input files (`text/`) |

---

## Installation

GRIFO needs Python with numpy, pandas and scikit-learn, Jupyter, and these command-line tools on
your `PATH`: NanoFilt, NanoStat, MultiQC, minimap2, spoa, vsearch, BLAST+, MAFFT and Krona. The easiest
way is one conda environment (Miniforge or Mambaforge):

```bash
mamba create -n grifo -c conda-forge -c bioconda \
    python=3.11 numpy pandas scikit-learn jupyterlab \
    nanofilt nanostat multiqc minimap2 spoa vsearch blast mafft krona
mamba activate grifo
```

Or create the same environment from the file in this folder:

```bash
mamba env create -f environment.yml
mamba activate grifo
```

Check that everything is found:

```bash
for t in NanoFilt NanoStat multiqc minimap2 spoa vsearch blastn makeblastdb mafft ktImportText; do
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

### 2. Edit `config.ini`

All settings are in `config.ini`, next to the notebook. Each one is explained in the file, with
the values allowed, the default, a suggestion and the effect of changing it. Before a first run,
check at least these:

| section | setting | what to set |
|---|---|---|
| `[paths]` | `path_data`, `path_to_blastdb` | your reads folder and the reference FASTA of your marker |
| `[NanoFilt]` | `minlength`, `maxlength`, `qscore` | read length range around your amplicon, and read quality |
| `[primers]` | `forward`, `reverse` | your primers, 5'→3' as ordered (IUPAC codes allowed) |
| `[clustering]` | `parallel_samples`, `threads` | together about the number of CPU cores |
| `[vsearch]`, `[haplonet]` | `threads` | the number of CPU cores |

After changing a setting, rerun the notebook section that uses it (named in each section heading
of `config.ini`) and the sections after it.

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
1. Every cluster sequence is turned to forward orientation (the strand of the forward primer),
   and the primers (with anything outside them) are cut off. If a primer occurs twice at an end
   (a primer concatemer), the inner copy is cut. Sequences without any primer are turned to forward
   orientation with `vsearch --orient`, using the others as reference. A primer is accepted with up to `max_mismatches` (4)
   mismatches, inserted or missing bases. Clusters with less than 30 bp between the primers (primer
   dimers) are left out, and with `require_both = yes` also clusters in which a primer is missing.
   The cell reports how many clusters and reads had both primers, one or none; a missing forward
   primer often means off-target amplification.
2. Identical sequences from all samples are merged into one ASV (amplicon sequence variant), with
   their reads added up. The ASVs are then clustered with vsearch into OTUs at `otu_identity`.
   The default, 1.0, only joins identical sequences, so each OTU is one ASV; set 0.97 for classic
   97% OTUs that group the variants of one species. OTUs are numbered by total reads (`OTU_1` is the
   most abundant).
3. Each OTU is BLASTed against the reference (`haplonet.py`). If the best hit leaves ≥ 10 bases at
   an end unmatched (a primer concatemer or chimeric piece), the OTU is trimmed to the matched part.
4. The taxonomy is the deepest rank on which ≥ 60% of the best hits agree. Ranks the identity
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
| `primers` | primers found in the OTU's sequence: `both`, `forward only`, `reverse only` or `none` |
| `n_asvs` | number of ASVs in the OTU |
| `sequence` | the OTU sequence (its most abundant ASV), between the primers, in forward orientation |

`asv_table.csv` has one row per ASV: its OTU, reads per sample, the OTU's taxonomy, the primers
found and the sequence; `asvs.fasta` has the ASV sequences.

`primer_trimming.csv` lists every sample cluster with the primers found, its orientation and its
length before and after trimming.

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
- **Min. best-hit identity** hides OTUs whose best BLAST hit is below the value you enter (in %),
  e.g. `98` to keep only OTUs with a close reference. OTUs without any hit are hidden for every
  value above 0; set it back to `0` to show all OTUs again. It combines with the text filter, and
  the line under it shows how many OTUs are left (e.g. *2 of 19 (filtered from 124, best hit ≥ 95%)*).
- The **red node** is your OTU. The other nodes are reference haplotypes, coloured by taxon;
  node size reflects the number of sequences, and the numbers on the edges are mutational steps.
- **Click a node** to list its accessions, identity to the OTU and full lineages.
- **Find** highlights the nodes matching a taxon or accession.
- **Below the network**, an NCBI-BLAST-style alignment shows the OTU against its best hit (the
  hit with the highest bit score, as on NCBI): score,
  E-value, identities, gaps, strand, query cover, and the alignment in blocks of 60 bases with
  mismatches and gaps in red. Pick another of the top hits from the list above it, or click a
  node to see the alignment of that haplotype. For MIDORI references the range is also given in
  the coordinates of the GenBank record, so you can compare it directly with NCBI web BLAST.
  Drag the line above the pane to resize it, or hide it.
- **Links for the best hit.** If the best hit is a eukaryote identified to species, the side panel
  and the alignment pane link to its **WoRMS** page, its **OBIS** distribution map and its **GBIF**
  occurrence map. The WoRMS and GBIF IDs are looked up online once per species, while section 6
  runs, and kept in `results/haplotype_networks/taxon_links.json`. Without internet, or for names
  WoRMS or GBIF don't know, the links search for the name instead (marked *search*; OBIS needs a
  WoRMS ID, so it is left out). If WoRMS lists the name as a synonym, the accepted name is shown.
- **Relax / Freeze layout** and dragging nodes tidy up the picture; **Save PNG** exports it.
- **Manual annotation** lets you replace GRIFO's assignment of an OTU, see below.

The side panel shows the OTU's read count, the OTU-table taxonomy and the consensus of the best
hits, plus a **Note** when the OTU was trimmed, had hits on N-runs discarded, or had no hits.

**Changing an assignment by hand.** When the network or the alignment shows that GRIFO's
assignment is wrong or too precise (for example, the OTU sits among the haplotypes of another
species than its best hit), correct it in the **Manual annotation** box of the side panel:

1. Under *Original* you see GRIFO's assignment. Fill in `scientificName` (the list suggests the
   taxa among the hits, but any name works) and `taxonRank`, and check `identifiedBy` and
   `dateIdentified` (your name is filled in from your previous annotation in this session, the
   date defaults to today).
2. Write why you changed it in `identificationRemarks`; it goes into the archive.
3. Click **Save for this OTU**. Changed OTUs are marked ✎ in the dropdown and show their new name;
   **Remove** takes the annotation away again. The OTU filter also finds manual names and remarks.
4. When you're done, click **Export manual_annotations.csv** and save the file as
   `metadata/manual_annotations.csv` (or where `[dwca] manual_annotations` points). Your browser
   saves it in your downloads folder; move it there yourself.

The export contains every annotation on the page: the ones from the existing
`manual_annotations.csv` (the page starts from that file) plus your edits, minus the ones you
removed.

The page doesn't remember anything itself: every time you open it, it starts from the
`manual_annotations.csv` it was built with, so edits from an earlier session can't slip into a
new file. Edits you haven't exported are lost when you close or reload the page (the browser asks
before it closes a page with unexported edits). Rerun section 6 after saving a new file, so that
the page starts from it. The
annotations change only the archive (section 7); the OTU table, the network page's own
assignment and the Krona plot keep GRIFO's assignment.

**Why lowering the identity threshold may not bring OTUs back.** The threshold can only show
OTUs that have a BLAST hit. With the default `blast_task = megablast`, BLAST finds only close
matches (above about 90% identity), so OTUs with only distant references get *no hits* at all,
even though `min_id` is 0.70. If many OTUs have no hits and every hit is at 98–100%, any threshold
between 1 and the lowest best-hit identity shows the same OTUs. To find more distant hits, set
`blast_task = dc-megablast` (or `blastn`, slower) in `[haplonet]` and run sections 5 and 6 again.
OTUs that still have no hits may be chimeras or off-target sequences (much longer or shorter than
your amplicon) or come from taxa missing from the reference.

**What the network tells you.** An OTU that sits on (0 steps) or right next to the haplotypes of
one species is a confident identification. An OTU between two species, or far from all of them,
points to a species missing from the reference, a sequencing artefact, or closely related
species the marker can't separate.

#### Section 7 – Darwin Core Archive
Packs the OTU table and your sample metadata into a Darwin Core Archive (DwC-A) for publishing
through a GBIF or OBIS IPT. The conversion is done by `edna2dwca.py` from
[MarGenoDaT](https://github.com/pascalhabluetzel/MarGenoDaT); the section prepares its input
from the GRIFO results.

**1. Make the sample metadata.** Open MarGenoDaT's `eDNA_metadata_template_builder.html` in a
browser, pick your fields (at least the *FAIRe mandatory* and *OBIS required* quick-add sets),
download `metadata_template.csv` and fill in one row per sample:

- `samp_name` must be the sample names of the OTU table, i.e. the folder names in `data/`
  (`barcode01`, …). Samples without a row are left out, with an error in the QC report.
- Mark blanks with `samp_category = negative control` (and `neg_cont_type`, e.g.
  `field negative`). Controls are left out of the occurrences.
- Give the primer sequences 5'→3' in `pcr_primer_forward` / `pcr_primer_reverse`; if they are
  missing, the sequences from `[primers]` are used.
- Comma, semicolon or tab separators and the builder's `#` guidance rows are all fine.

`metadata/sample_metadata_mock.csv` shows a filled-in example. It is **made up** (stations,
dates and sampling details are invented) and only there to test the section.

**2. Fill in `[dwca]` in `config.ini`:** the path to the metadata, and the dataset title,
abstract, licence and creators for the EML. Only the first creator is required. Four optional
settings hold what GRIFO can't derive from the data: `target_subfragment`, `pcr_primer_reference`
(GBIF and OBIS recommend it; the QC report warns while it is empty), `adapters` and
`associated_sequences` (the ENA/SRA accession of the raw reads, once deposited). Empty values are
left out of the archive.

**3. Run the section.** It writes to `results/dwca/`:

| file | contents |
|---|---|
| `GRIFO_dwca.zip` | the archive: `occurrence.txt`, `dna_derived_data.txt`, `meta.xml`, `eml.xml` |
| `GRIFO_dwca_QC_report.md` | errors, warnings and inferred values; also shown in the notebook |
| `otu_table_dwca.tsv`, `otus_dwca.fasta` | the OTU table and sequences as `edna2dwca.py` reads them |
| `dataset_config.json` | the EML information, pipeline defaults and per-sample values passed to `edna2dwca.py` |

The section also writes `results/qc/raw_read_counts.csv`: the raw reads per sample, counted once
from the FASTQ files in `data/` (this takes a moment the first time) and reused afterwards.

Each record is one OTU in one sample with reads: `organismQuantity` is its reads,
`sampleSizeValue` the sample's total reads. The taxonomy is the OTU table's (after the identity
cutoffs). The lineage, with NCBI taxon IDs as `taxonID`, comes from the MIDORI2 header of the best
hit, or of the first top hit that belongs to the consensus taxon. `verbatimIdentification` is
GRIFO's call (e.g. `Chelon auratus (species)`), and `scientificNameID` the WoRMS LSID of the name
when section 6 found one (`[haplonet] taxon_links = yes`). Unassigned OTUs are kept as *Biota*
with their sequence (`[dwca] unassigned = drop` leaves them out).

**Manual identifications.** OTUs in `metadata/manual_annotations.csv` (see section 6; the file
is optional) are published with your identification, and GRIFO's is kept:

| term | value |
|---|---|
| `scientificName`, `taxonRank`, lineage, `taxonID` | your name and rank; lineage and NCBI taxon ID from the first MIDORI2 reference whose lineage contains the name |
| `scientificNameID` | the WoRMS LSID of your name (looked up online when `[haplonet] taxon_links = yes`) |
| `identifiedBy`, `dateIdentified` | from the file |
| `identificationVerificationStatus` | `manually verified` |
| `identificationRemarks` | your remarks, followed by *Manually identified by …; automatic assignment in previousIdentifications* |
| `previousIdentifications` | GRIFO's assignment, e.g. `Chelon auratus (species), 100.0% identity of the best hit (BLAST+ megablast best-hit consensus, GRIFO v0.8)` |
| `verbatimIdentification` | GRIFO's call, unchanged (`no hits` for unassigned OTUs) |
| FAIRe `percent_match`, `confidence_score`, … | unchanged: they describe the BLAST result, not your identification |

A name that isn't in the reference lineages (a typo, or a taxon missing from MIDORI2) is
published without taxon ID and higher ranks; the notebook and the QC report warn about it.
Annotations of OTUs that aren't in the OTU table are ignored with a warning. The QC report lists
all manual identifications with their previous assignment.

Everything GRIFO knows about the analysis goes into the archive under an official term. Where
the sample metadata has a value for a term, the metadata wins.

| what | where it comes from | terms |
|---|---|---|
| sample and assay fields | sample metadata | every Darwin Core and DNA-derived data extension column of the metadata |
| OTU sequence, primers, amplicon size | `otus.fasta`, `[primers]`, OTU length + primers (if both were found) | `DNA_sequence`, `pcr_primer_forward/reverse`, `ampliconSize` |
| read counts per sample | FASTQ files, OTU table | `lib_reads_seqd`; FAIRe `input_read_count`, `output_read_count`, `output_otu_num`, `otu_num_tax_assigned` |
| barcode | sample name (`barcode01` …) | `mid` (the native barcode name, not its sequence) |
| read QC | `[NanoFilt]` | `seq_quality_check`; FAIRe `min_len_cutoff`, `min_len_tool`, `error_rate_tool/type/cutoff` |
| primer trimming | `[primers]` | FAIRe `trim_method`, `trim_param` |
| clustering and read assignment | `[clustering]`, `[read_assignment]`, `[vsearch]` | `otu_class_appr`, `chimera_check`; FAIRe `otu_clust_tool`, `otu_clust_cutoff`, `chimera_check_method`, `bioinfo_method_additional` |
| BLAST and taxonomy | `[haplonet]` | `otu_db`, `ref_db`, `otu_seq_comp_appr`, `sim_search_meth`, `tax_class`; FAIRe `tax_assign_cat`, `tax_class_id_cutoff`, `tax_class_query_cutoff`, `tax_class_collapse`, `tax_class_other`, `taxonID_db` |
| assignment of each OTU | OTU table | FAIRe `percent_match`, `percent_query_cover`, `confidence_score`, `accession_id_ref_db`, `seq_id` |
| filtering and screening | `[dwca] min_reads`, QC report | FAIRe `min_reads_cutoff`, `min_reads_cutoff_unit`, `min_reads_tool`, `screen_nontarget_method`, `screen_contam_method` |
| pipeline | `[dwca] identification_reference` | `sop`, `identificationReferences`; FAIRe `code_repo`, `sop_bioinformatics` |

Darwin Core and DNA-derived data extension terms get their own columns. FAIRe terms that are
neither (the FAIRe checklist isn't part of the extension yet) go into `dynamicProperties` as
JSON under their FAIRe names, and the dataset-wide ones are also listed in the EML methods. Tool
versions are read from the programs on your PATH when the section runs.

Fix every **error** in the QC report and read the warnings (for example OTUs outside the target
taxon, or a missing primer reference). Then check the archive with the
[GBIF data validator](https://www.gbif.org/tools/data-validator) before uploading it to an IPT.

`edna2dwca.py` is MarGenoDaT's script with some changes, listed at the top of the file:
semicolon files and `#` guidance rows are read, the FAIRe field `samp_collect_method` is used,
the identity is described as the best-hit identity, all metabarcoding terms of the DNA-derived
data extension and the FAIRe terms (in `dynamicProperties`, with the OTU ID as `seq_id`) are
written, the OTU table may carry per-OTU term columns (including manual identifications), and
the config file may hold per-sample values (`sample_values`). It still runs on its own: `python edna2dwca.py --help`.

#### Section 8 – Krona plot
Draws the OTU-table taxonomy as an interactive [Krona](https://github.com/marbl/Krona) chart,
`results/krona/krona.html`: wedges from domain to species, sized by reads, with the OTUs as the
outermost ring, named with their taxon and best-hit identity (*Chelon auratus - OTU_1 (100.0%)*).
Click a wedge to zoom in, or hover over it for its name; the chart works offline. Krona's
*Collapse* box starts unticked (`[krona] collapse = no`), so every rank has its own ring.
Inner rings get thin when there are many ranks: zoom into a wedge to read genus and species, or
set `otu_level = no` to make the deepest taxon the outermost ring.

**Identity threshold.** A Krona chart can't filter itself, so it holds one dataset per value in
`[krona] thresholds` (default 0, 80, 85, 90, 95, 97, 98, 99, 100). Pick one in the list at the top
left, or step through them with ↑ ↓. Each leaves out the OTUs whose best BLAST hit is below the
threshold, as *Min. best-hit identity* does on the network page. OTUs without hits appear as
*No hit* at 0% only. The dataset names give the OTUs and reads left, e.g.
*All samples, min. identity 98% (19 OTUs, 180,310 reads)*. With `per_sample = yes` (default)
there is the same set for each sample; control samples are included.

Needs `ktImportText` from the conda package `krona`. Only `ktImportText` is used, so you
don't need to download Krona's NCBI taxonomy (`ktUpdateTaxonomy.sh`), which conda suggests after
the install. To add Krona to an existing environment:

```bash
mamba install -n grifo -c conda-forge -c bioconda krona
```

Options in `[krona]` (`config.ini`):

```ini
thresholds = 0, 80, 85, 90, 95, 97, 98, 99, 100   ; minimum best-hit identity (%), one dataset each
per_sample = yes      ; also one set of datasets per sample (no = all samples together only)
otu_level = yes       ; OTUs as the outermost ring (no = the deepest taxon, OTUs added up)
collapse = no         ; start with Krona's 'Collapse' ticked (merges ranks with a single child)
```

To make the chart again after changing them, run the first code cell of section 5 and section 8.

### 4. Check a result at NCBI

To compare an OTU with NCBI's web BLAST, copy its `sequence` from `otu_table.csv` (or `otus.fasta`)
of the **same run**. OTU numbers are assigned by read count in each run, so `OTU_1` today can be a
different sequence from `OTU_1` in an earlier run or results folder.

Keep in mind that NCBI searches all of GenBank, while GRIFO only searches your reference file. An OTU
with *no hits* here may be off-target DNA (bacteria, nuclear DNA) that a mitochondrial database
cannot identify.

---

## Settings

Every setting is in `config.ini`, each with an explanation, the values allowed, the default, a
suggestion and the effect of changing it. Sections:

| section | used by | settings |
|---|---|---|
| `[paths]` | all | reads folder, results folder, reference FASTA |
| `[NanoFilt]` | 2 | read length range, read quality |
| `[clustering]` | 4 | iterations, parallel samples and threads, consensus and clustering parameters, random seed |
| `[read_assignment]` | 4 | identity and coverage a read needs to be counted, minimap2 preset |
| `[primers]` | 5 | primers, allowed mismatches, removing sequences without both primers |
| `[vsearch]` | 5 | identity for clustering ASVs into OTUs, threads |
| `[haplonet]` | 5, 6 | BLAST database and search, taxonomy cutoffs, haplotype networks and page |
| `[dwca]` | 6, 7 | sample metadata file, manual annotations file, dataset description (EML), unassigned OTUs, minimum reads, primer reference, target subfragment, adapters, raw-read accession |
| `[krona]` | 8 | identity thresholds of the datasets, per-sample charts, OTU ring |
| `[BLAST]` | – | only used by the old `grifo.py` script |

If you delete a line, its default is used (except the required settings in `[paths]`,
`[NanoFilt]` and `[clustering] niter`).

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
`hn_*` columns added. Add `--annotations metadata/manual_annotations.csv` to start the page's
**Manual annotation** box from an existing annotation file.

---

## Troubleshooting

| problem | likely cause / fix |
|---|---|
| `Required tool(s) not found on PATH` | the conda environment isn't active in the Jupyter kernel; start Jupyter from the activated environment |
| a sample called `data` appears | FASTQ files at the top level of `data/`; move them into sample folders |
| reads counted twice | merged FASTQ files next to the original chunks in a sample folder |
| empty ranks in the header check (section 5) | reference header format not recognised; use a supported format or a `--taxonomy` TSV with `haplonet.py` |
| `BLAST hit IDs … are not in …` | `blastdb` was built from a different FASTA than `db` |
| few clusters with `both` primers | wrong primer sequences, or the reverse primer entered reverse-complemented; enter both 5'→3' as ordered |
| many `unassigned` OTUs | wrong marker database, off-target amplification, or too strict cutoffs; look at `identity` and `consensus`, try `blast_task = dc-megablast` |
| lowering *Min. best-hit identity* on the network page doesn't bring OTUs back | those OTUs have no BLAST hit at all; megablast finds only close matches, so set `blast_task = dc-megablast` and run sections 5 and 6 again |
| species missing although `consensus` names one | the best identity is below `species_identity` |
| network page is blank | no internet connection (vis-network loads from a CDN) |
| `Sample '…' is in the OTU table but has no metadata row` (section 7) | `samp_name` in the metadata doesn't match the sample (folder) names; use `barcode01`, … |
| `Required tool(s) not found on PATH: ktImportText` (section 8) | Krona isn't installed: `mamba install -n grifo -c bioconda krona` |
| manual annotations gone from the network page | the page was closed or reloaded before you exported them, or it was built (section 6) before you saved the new `manual_annotations.csv`; rerun section 6 after saving it |
| `Not found in the reference lineages` / `Manual identifications without NCBI taxon ID` (section 7) | the manual name has a typo, or another rank in MIDORI2 than the `taxonRank` you gave (leave `taxonRank` empty to accept any rank), or the taxon is missing from the reference |
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
