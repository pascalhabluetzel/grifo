#!/usr/bin/env python3
"""
haplonet.py - Taxonomic assignment + haplotype networks for many sequences
(e.g. every OTU of an OTU table), shown in ONE interactive HTML page.

For every sequence (OTU):
  1. BLAST (blastn, default) or VSEARCH (--search vsearch) searches it
     against a reference database and reports a taxonomic assignment
     (consensus of the best hits).
  2. The best N hits (default 100) are kept. If the best BLAST hit leaves
     >= 10 bases at an OTU end unmatched (e.g. a primer concatemer or chimeric
     piece), the OTU is trimmed to the matched part. Hits covering < 80% of
     the (trimmed) OTU and hits on N-runs in the reference are discarded.
  3. The OTU + its hits are aligned with MAFFT (only the region of each
     reference that matches the OTU).
  4. The alignment is trimmed to the span covered by the OTU.
  5. Identical sequences are collapsed into haplotypes and a minimum spanning
     network is built. Haplotypes are labelled with their taxonomy (species if
     available, otherwise the lowest available rank).

All networks go into one HTML file: pick an OTU from the dropdown menu, or
step through the OTUs with the arrow keys (left/right or up/down). Below each
network, an NCBI-BLAST-style alignment shows the OTU against its best hit
(--search blast); the other top hits can be picked from a list, or by clicking
a node in the network.

Input (-i / --input)
  - an OTU table (CSV or TSV) with an OTU-ID column and a sequence column
    (auto-detected: "OTU"/"id"/first column, and "sequence"/"seq"), or
  - a FASTA file with one record per OTU, or
  - a single raw DNA sequence.
  Other columns of a table (read counts, an existing taxonomy) are shown in
  the page and kept in the output table.

Requirements
  - Python >= 3.8 (standard library only)
  - BLAST+ (blastn, makeblastdb) and mafft on the PATH
      conda install -c bioconda blast mafft
    (vsearch only if you use --search vsearch)
  - The HTML file loads vis-network from a CDN (internet needed to view it).
  - For best hits that are eukaryote species, the page links to WoRMS, OBIS
    and GBIF. The WoRMS and GBIF IDs are looked up online once per species
    (cached in <outdir>/taxon_links.json); offline, or with
    --no-taxon-links, the links search for the name instead.

Supported reference header formats (taxonomy parsed automatically)
  MIDORI2       : >ID root_1;...;genus_Homo_9605;species_Homo_sapiens_9606
                  (also with '###' between ID and lineage)
  SINTAX / UTAX : >ID;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_subtilis;
  QIIME / UNITE : >ID k__Fungi;p__Ascomycota;...;s__Penicillium_roqueforti
  SILVA / PR2   : >ID Bacteria;Firmicutes;Bacilli;...;Bacillus subtilis
  Plain         : >ID Genus species
  Or give a separate TSV (ID <tab> taxonomy string) with --taxonomy.

BLAST database
  -d is always the reference FASTA (sequences and taxonomy are read from it).
  Without --blastdb the BLAST database is built from it in <outdir>/blastdb;
  reuse that with --blastdb to skip the build next time. To build it yourself:
      sed 's/###/ /' MIDORI2.fasta > MIDORI2_blast.fasta     # MIDORI headers only
      makeblastdb -in MIDORI2_blast.fasta -dbtype nucl -out MIDORI2_lrRNA
  (do not use -parse_seqids: MIDORI IDs are not valid BLAST seq-ids)

Examples
  python haplonet.py -i otu_table.csv -d MIDORI2_lrRNA.fasta -o results
  python haplonet.py -i otu_table.csv -d MIDORI2_lrRNA.fasta -o results2 \\
         --blastdb results/blastdb/reference
  python haplonet.py -i otu_table.csv -d ref.fasta -o results --search vsearch
  python haplonet.py -i otus.fasta    -d reference.fasta -o results --threads 8
  python haplonet.py -i otu_table.csv -d ref.fasta -o results --otus OTU_1,OTU_3
"""

import argparse
import concurrent.futures as cf
import csv
import html
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, OrderedDict, defaultdict

# --------------------------------------------------------------------------
# FASTA helpers
# --------------------------------------------------------------------------

def read_fasta(path):
    """Return an OrderedDict {full_header: sequence}."""
    seqs = OrderedDict()
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        seqs[header] = "".join(chunks)
    return seqs


def write_fasta(records, path, width=80):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


_COMP = str.maketrans("ACGTURYKMBDHVNacgturykmbdhvn-",
                      "TGCAAYRMKVHDBNtgcaayrmkvhdbn-")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def seq_id(header):
    """Sequence ID: header up to the first whitespace, ';' or '|'."""
    return re.split(r"###|[\s;|]", header.strip(), maxsplit=1)[0] or header


# --------------------------------------------------------------------------
# Taxonomy parsing
# --------------------------------------------------------------------------

RANK_ORDER = ["domain", "kingdom", "phylum", "class", "order",
              "family", "genus", "species"]
RANK_CODES = {"d": "domain", "k": "kingdom", "p": "phylum", "c": "class",
              "o": "order", "f": "family", "g": "genus", "s": "species"}

RANK_WORDS = {"superkingdom": "domain", "realm": "domain", "domain": "domain",
              "kingdom": "kingdom", "phylum": "phylum", "class": "class",
              "order": "order", "family": "family", "genus": "genus",
              "species": "species"}
MIDORI_RANK = re.compile(
    r"^(root|cellularroot|superkingdom|realm|domain|kingdom|subkingdom|"
    r"superphylum|phylum|subphylum|superclass|class|subclass|infraclass|"
    r"superorder|order|suborder|infraorder|superfamily|family|subfamily|"
    r"tribe|subtribe|genus|subgenus|species.group|species.subgroup|species|clade|cohort|subcohort|section|"
    r"subsection|series|parvorder|infraphylum|subspecies|varietas|forma|"
    r"no.rank)_(.+)$", re.IGNORECASE)

UNINFORMATIVE = re.compile(
    r"^(|na|n/a|none|null|unknown|unidentified|unclassified|uncultured|"
    r"incertae.?sedis|metagenome|environmental.?sample|.*_sp\.?|.* sp\.?|"
    r"uncultured.*|unidentified.*|unclassified.*)$",
    re.IGNORECASE,
)


def clean_name(name):
    name = name.strip().strip('"').strip("'").replace("_", " ")
    name = re.sub(r"\s+", " ", name).strip()
    name = re.sub(r"\s\d+$", "", name)          # trailing taxid (MIDORI)
    return name


def informative(name):
    return bool(name) and not UNINFORMATIVE.match(name.strip())


def parse_taxonomy_string(tax):
    """
    Parse a taxonomy string into a list of (rank, name) from high to low.
    rank may be None when the format carries no rank prefixes.
    """
    tax = tax.strip().rstrip(";").strip()
    if not tax:
        return []

    # SINTAX: d:Bacteria,p:Firmicutes,...
    if re.match(r"^[a-z]:", tax) and ("," in tax or tax.count(":") == 1):
        out = []
        for part in tax.split(","):
            if ":" in part:
                code, name = part.split(":", 1)
                out.append((RANK_CODES.get(code.strip().lower()), clean_name(name)))
        return out

    parts = [p for p in re.split(r"[;|]", tax) if p.strip()]
    out = []
    for p in parts:
        p = p.strip()
        m = re.match(r"^([a-z])__(.*)$", p)            # QIIME / UNITE / GTDB
        if m:
            out.append((RANK_CODES.get(m.group(1)), clean_name(m.group(2))))
            continue
        m = MIDORI_RANK.match(p)                         # MIDORI: species_Homo_sapiens_9606
        if m:
            word = m.group(1).lower()
            if word in ("root", "cellularroot", "noRank".lower()):
                continue
            out.append((RANK_WORDS.get(word, word), clean_name(m.group(2))))
            continue
        m = re.match(r"^([a-z]):(.*)$", p)             # rank:name without commas
        if m:
            out.append((RANK_CODES.get(m.group(1)), clean_name(m.group(2))))
            continue
        out.append((None, clean_name(p)))
    return out


def taxonomy_from_header(header):
    """Extract a taxonomy string from a reference header (best effort)."""
    m = re.search(r"tax=([^;]+(?:;[^;=]+)*?);?\s*$", header)
    if m:
        return m.group(1)
    m = re.search(r"tax=([^;]+)", header)
    if m:
        return m.group(1)

    if "###" in header:
        return header.split("###", 1)[1]
    parts = header.split(None, 1)
    rest = parts[1] if len(parts) > 1 else ""
    if not rest and ";" in header:              # ID;lineage;...
        rest = header.split(";", 1)[1]

    # BOLD-like: ID|Genus species|marker|...
    if not rest and "|" in header:
        fields = header.split("|")[1:]
        for f in fields:
            if re.match(r"^[A-Z][a-z]+([ _][a-z\-]+)?$", f.strip()):
                return f
        return ";".join(fields)
    return rest


def lowest_label(lineage):
    """
    Species if available (binomial if the species name is only an epithet),
    otherwise the lowest informative rank. Returns (label, rank).
    """
    if not lineage:
        return "unassigned", None
    ranked = {r: n for r, n in lineage if r}
    sp = ranked.get("species")
    if sp and informative(sp):
        genus = ranked.get("genus")
        if genus and " " not in sp and not sp.startswith(genus):
            sp = f"{genus} {sp}"
        return sp, "species"
    for rank, name in reversed(lineage):
        if rank == "species" or (rank and rank not in RANK_ORDER):
            continue
        if informative(name):
            if rank is None and re.match(r"^[A-Z][a-z]+ [a-z][a-z\-]+", name):
                return name, "species"
            return name, rank or "lowest"
    return "unassigned", None


def load_taxonomy_tsv(path):
    tax = {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and not row[0].startswith("#"):
                tax[row[0].strip()] = row[1].strip()
    return tax


# --------------------------------------------------------------------------
# External tools
# --------------------------------------------------------------------------

def run(cmd, log=None):
    print("  $ " + " ".join(cmd), file=sys.stderr)
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         text=True)
    if log:
        with open(log, "a") as fh:
            fh.write("$ " + " ".join(cmd) + "\n" + res.stderr + "\n")
    if res.returncode != 0:
        sys.exit(f"ERROR running {cmd[0]}:\n{res.stderr}")
    return res


def check_tool(name, path):
    exe = shutil.which(path)
    if not exe:
        sys.exit(f"ERROR: '{path}' not found. Install {name} "
                 f"(e.g. conda install -c bioconda {name}) or pass its path.")
    return exe


# --------------------------------------------------------------------------
# Step 1-2: VSEARCH search and taxonomic assignment
# --------------------------------------------------------------------------


def vsearch_hits(vsearch, query_fa, db_fa, outdir, top, min_id, threads,
                 strand, maxrejects, log):
    """
    Search all queries at once; return {query_label: [hits sorted best first]}.
    `top` here is the number of hits requested from VSEARCH (already including
    spare hits, because hits on N-runs are discarded afterwards).
    """
    out_tsv = os.path.join(outdir, "vsearch_hits_raw.tsv")
    fields = "query+target+id+alnlen+mism+opens+qcov+tcov+ql+tl+qstrand+tilo+tihi"
    cmd = [vsearch, "--usearch_global", query_fa, "--db", db_fa,
           "--id", str(min_id), "--maxaccepts", str(top),
           "--maxrejects", str(maxrejects), "--maxhits", str(top),
           "--strand", strand, "--userout", out_tsv, "--userfields", fields,
           "--threads", str(threads), "--notrunclabels"]
    run(cmd, log)
    names = fields.split("+")
    by_query = defaultdict(list)
    with open(out_tsv) as fh:
        for line in fh:
            row = dict(zip(names, line.rstrip("\n").split("\t")))
            row["id"] = float(row["id"])
            by_query[row["query"]].append(row)
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(names) + "\n")
        for q in by_query:
            hits = sorted(by_query[q], key=lambda h: (-h["id"], -int(h["alnlen"])))
            by_query[q] = hits
            for h in hits:
                fh.write("\t".join(str(h[n]) for n in names) + "\n")
    return by_query


def make_blast_db(makeblastdb, db_fa, outdir, log):
    """
    Build a BLAST nucleotide database from the reference FASTA. MIDORI headers
    (ID###lineage) are rewritten to 'ID lineage' so BLAST sees a short ID.
    -parse_seqids is deliberately not used: MIDORI IDs such as
    'MT217036.1.<1.>530' are not valid BLAST seq-ids.
    """
    dbdir = os.path.join(outdir, "blastdb")
    os.makedirs(dbdir, exist_ok=True)
    prefix = os.path.join(dbdir, "reference")
    src = db_fa
    with open(db_fa) as fh:
        needs_fix = any("###" in line for _, line in zip(range(2000), fh)
                        if line.startswith(">"))
    if needs_fix:
        src = os.path.join(dbdir, "reference_clean.fasta")
        with open(db_fa) as fi, open(src, "w") as fo:
            for line in fi:
                fo.write(line.replace("###", " ", 1) if line.startswith(">") else line)
    run([makeblastdb, "-in", src, "-dbtype", "nucl", "-out", prefix,
         "-title", os.path.basename(db_fa)], log)
    if needs_fix:
        os.remove(src)
    print(f"      BLAST database written to {prefix} "
          f"(reuse it next time with --blastdb {prefix})", file=sys.stderr)
    return prefix


def blast_hits(blastn, query_fa, blastdb, outdir, top, min_id, threads, task,
               evalue, min_qcov, log):
    """
    blastn all queries at once; return {query_label: [hits sorted best first]}
    with the same keys as vsearch_hits(). Each hit is the best HSP per
    reference. tilo/tihi are extended over the query ends BLAST left out
    (local alignment) later, in fit_hits_to_query().
    """
    out_tsv = os.path.join(outdir, "blast_hits_raw.tsv")
    cols = ("qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
            "evalue bitscore score qlen slen qcovhsp qseq sseq")
    cmd = [blastn, "-task", task, "-query", query_fa, "-db", blastdb,
           "-max_target_seqs", str(top), "-max_hsps", "1",
           "-perc_identity", str(round(min_id * 100, 2)), "-evalue", str(evalue),
           "-num_threads", str(threads), "-outfmt", "6 " + cols, "-out", out_tsv]
    run(cmd, log)
    names = cols.split()
    by_query = defaultdict(dict)
    with open(out_tsv) as fh:
        for line in fh:
            r = dict(zip(names, line.rstrip("\n").split("\t")))
            qs, qe, ss, se = (int(r[k]) for k in ("qstart", "qend", "sstart", "send"))
            ql, sl = int(r["qlen"]), int(r["slen"])
            plus = ss <= se
            h = {"query": r["qseqid"], "target": r["sseqid"],
                 "id": float(r["pident"]), "alnlen": int(r["length"]),
                 "mism": int(r["mismatch"]), "opens": int(r["gapopen"]),
                 "qcov": float(r["qcovhsp"]),
                 "tcov": round(100 * int(r["length"]) / sl, 1),
                 "ql": ql, "tl": sl, "qstrand": "+" if plus else "-",
                 "qstart": qs, "qend": qe, "sstart": ss, "send": se,
                 "evalue": r["evalue"], "bitscore": float(r["bitscore"]),
                 "score": int(float(r["score"])), "qaln": r["qseq"], "saln": r["sseq"]}
            old = by_query[h["query"]].get(h["target"])
            if old is None or h["bitscore"] > old["bitscore"]:
                by_query[h["query"]][h["target"]] = h
    return {q: sorted(d.values(), key=lambda h: (-h["id"], -h["qcov"], -h["bitscore"]))
            for q, d in by_query.items()}


def fit_hits_to_query(hits, qlen, min_qcov, trim, trim_min):
    """
    BLAST aligns locally. Find the part of the query the best hit(s) cover
    (the 'core'); if `trim` and more than `trim_min` bases lie outside it,
    the query is trimmed to the core (e.g. a primer concatemer or chimeric
    piece at one end). Query cover is then measured against the core, hits
    below `min_qcov` are dropped, and each reference region (tilo..tihi) is
    extended/clipped so that it spans exactly the core.
    Returns (hits, core_lo, core_hi) with 1-based inclusive core coordinates.
    """
    lo, hi = 1, qlen
    if hits and trim:
        best = max(h["bitscore"] for h in hits)
        top = [h for h in hits if h["bitscore"] >= best]
        clo, chi = min(h["qstart"] for h in top), max(h["qend"] for h in top)
        if (clo - 1) + (qlen - chi) >= trim_min:
            lo, hi = clo, chi
    clen = hi - lo + 1
    kept = []
    for h in hits:
        qs, qe, ss, se, sl = h["qstart"], h["qend"], h["sstart"], h["send"], h["tl"]
        overlap = max(0, min(qe, hi) - max(qs, lo) + 1)
        h["qcov"] = round(100 * overlap / clen, 1)
        if h["qcov"] < min_qcov:
            continue
        if ss <= se:
            tlo, thi = ss - (qs - lo), se + (hi - qe)
        else:
            tlo, thi = se - (hi - qe), ss + (qs - lo)
        h["tilo"], h["tihi"] = max(1, tlo), min(sl, thi)
        kept.append(h)
    kept.sort(key=lambda h: (-h["id"], -h["qcov"], -h["bitscore"]))
    return kept, lo, hi


def consensus_label(lineages, cutoff):
    """
    Rank-by-rank consensus: descend while the most common taxon at a rank is
    shared by >= `cutoff` of the hits (1.0 = strict lowest common ancestor).
    Returns (label, rank, support).
    """
    dicts = []
    for lin in lineages:
        d = {r: n for r, n in lin if r in RANK_ORDER and informative(n)}
        if "species" in d and "genus" in d and " " not in d["species"]:
            d["species"] = f"{d['genus']} {d['species']}"
        dicts.append(d)
    if not any(dicts):                              # unranked lineages
        dicts = [{i: n for i, (_, n) in enumerate(lin) if informative(n)}
                 for lin in lineages]
        ranks = sorted({k for d in dicts for k in d})
    else:
        ranks = RANK_ORDER
    label, lrank, support = "unassigned", None, 0.0
    N = len(dicts)
    for rank in ranks:
        names = [d.get(rank) for d in dicts]
        present = [x for x in names if x]
        if not present:
            continue                                # rank not used in this DB
        top, k = Counter(present).most_common(1)[0]
        if k / N < cutoff:                          # missing values count against
            break
        label, lrank, support = top, rank, k / N
    return label, (lrank if isinstance(lrank, str) else "lowest"), support


def consensus_assignment(hits, tax_of, delta=0.0, cutoff=0.6):
    """
    Take the hits at the best identity (or within `delta` %-points of it) and
    report the lowest taxon supported by >= `cutoff` of them.
    """
    if not hits:
        return None
    best = hits[0]["id"]
    tied = [h for h in hits if h["id"] >= best - delta]
    label, rank, support = consensus_label(
        [tax_of[h["target"]]["lineage"] for h in tied], cutoff)
    labels = Counter(tax_of[h["target"]]["label"] for h in tied)
    return {"label": label, "rank": rank, "support": support, "identity": best,
            "n_top_hits": len(tied), "top_targets": [h["target"] for h in tied],
            "top_labels": labels}


def alignment_records(hits, label_of, n=10):
    """
    Pairwise BLAST alignments of the OTU with its best `n` hits by bit score
    (the order NCBI uses), for the NCBI-style alignment pane of the HTML page. `label_of(target)` gives the
    taxon name shown. Coordinates are those BLAST reports: on the untrimmed OTU
    and on the reference sequence. Hits without a stored alignment (VSEARCH
    search) are skipped.
    """
    out = []
    ranked = sorted(hits, key=lambda h: -(h.get("bitscore") or 0))   # stable: ties keep their order
    for h in ranked[:max(0, n)]:
        q, s = (h.get("qaln") or "").upper(), (h.get("saln") or "").upper()
        if not q or len(q) != len(s):
            continue
        out.append({"target": h["target"], "acc": seq_id(h["target"]),
                    "label": label_of(h["target"]),
                    "qs": h["qstart"], "qe": h["qend"], "ql": h["ql"],
                    "ss": h["sstart"], "se": h["send"], "sl": h["tl"],
                    "bits": h.get("bitscore"), "score": h.get("score"),
                    "ev": h.get("evalue"), "pid": h["id"],
                    "ident": sum(a == b and a != "-" for a, b in zip(q, s)),
                    "len": len(q), "gaps": q.count("-") + s.count("-"),
                    "q": q, "s": s})
    return out


# --------------------------------------------------------------------------
# Links to WoRMS, OBIS and GBIF for the best hit
# --------------------------------------------------------------------------
EUKARYOTE_KINGDOMS = {"animalia", "metazoa", "plantae", "viridiplantae", "fungi", "chromista",
                      "protozoa", "protista", "stramenopiles", "alveolata", "rhizaria",
                      "amoebozoa", "excavata", "haptophyta", "cryptophyta", "rhodophyta",
                      "chlorophyta", "discoba", "metamonada", "sar"}
SPECIES_NAME = re.compile(r"^[A-Z][a-z]+ [a-z][a-z\-]+( [a-z][a-z\-]+)?$")
WORMS_URL = "https://www.marinespecies.org/rest/AphiaRecordsByName/{}?like=false&marine_only=false"
# GBIF switched its taxonomy to the Catalogue of Life eXtended Release (COL XR): taxon keys are
# now short codes (e.g. Q2M4); the old integer backbone keys return no occurrences on gbif.org
GBIF_COL_XR = "7ddf754f-d193-4cc9-b351-99906754a03b"
GBIF_URL = "https://api.gbif.org/v2/species/match?scientificName={}&checklistKey=" + GBIF_COL_XR
GBIF_COUNT_URL = ("https://api.gbif.org/v1/occurrence/search?taxonKey={}&checklistKey="
                  + GBIF_COL_XR + "&limit=0")
GBIF_MAP_URL = "https://www.gbif.org/occurrence/map?taxonKey={}"
WORMS_SYN_URL = "https://www.marinespecies.org/rest/AphiaSynonymsByAphiaID/{}?offset=1"


def best_hit(hits):
    """The hit with the highest bit score (the 'Best hit' of the HTML page), or None."""
    return max(hits, key=lambda h: h.get("bitscore") or 0) if hits else None


def link_species(tax):
    """Species name of a reference that is a eukaryote identified to species ({'lineage',
    'label', 'rank'} as from lowest_label), for the WoRMS/OBIS/GBIF links; else None."""
    if not tax or tax.get("rank") != "species":
        return None
    ranks = {r: n.lower() for r, n in tax["lineage"] if r}
    eukaryote = (ranks.get("domain", "").startswith("eukaryot")
                 or ranks.get("kingdom") in EUKARYOTE_KINGDOMS)
    return tax["label"] if eukaryote and SPECIES_NAME.match(tax["label"]) else None


def _get_json(url, timeout=15):
    import urllib.request
    req = urllib.request.Request(url, headers={"User-Agent": "haplonet (GRIFO pipeline)",
                                               "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        body = r.read()
    return json.loads(body) if body.strip() else None


def _gbif_key(name):
    """GBIF (COL XR) taxon key of the accepted species for an exact species-name match, else
    None. Reads the v2 match response ({usage, acceptedUsage, diagnostics}), and a flat one."""
    from urllib.parse import quote
    m = _get_json(GBIF_URL.format(quote(name))) or {}
    usage = m.get("usage") or {}
    match_type = (m.get("diagnostics") or {}).get("matchType") or m.get("matchType")
    rank = str(usage.get("rank") or m.get("rank") or "").upper()
    key = ((m.get("acceptedUsage") or {}).get("key") or usage.get("key")
           or m.get("acceptedUsageKey") or m.get("usageKey"))
    if match_type == "EXACT" and rank in ("SPECIES", "SUBSPECIES") and key:
        return str(key)
    return None


def _gbif_best_key(names, seen):
    """(key, occurrences) of the GBIF key with the most occurrences among the names."""
    best = (None, -1)
    for n in names:
        key = _gbif_key(n)
        if not key or key in seen:
            continue
        seen.add(key)
        count = int((_get_json(GBIF_COUNT_URL.format(key)) or {}).get("count") or 0)
        if count > best[1]:
            best = (key, count)
    return best


def _lookup_species(name):
    """Accepted WoRMS AphiaID (+ accepted name), and the GBIF species key with occurrences, of a
    species name; None where not found. GBIF may file a species' occurrences under another
    name than WoRMS accepts (e.g. Liza aurata for Chelon auratus): if neither the name nor the
    WoRMS accepted name has occurrences in GBIF, the WoRMS synonyms are tried, and the key with
    the most occurrences is used. Raises on network errors."""
    from urllib.parse import quote
    out = {"aphia": None, "aphia_name": None, "gbif": None, "gbif_n": None, "gbif_checklist": None}
    records = _get_json(WORMS_URL.format(quote(name))) or []
    exact = [r for r in records if (r.get("scientificname") or "").lower() == name.lower()]
    if exact:
        r = next((x for x in exact if x.get("status") == "accepted"), exact[0])
        out["aphia"] = r.get("valid_AphiaID") or r.get("AphiaID")
        out["aphia_name"] = r.get("valid_name") or r.get("scientificname")
    seen = set()
    key, count = _gbif_best_key(dict.fromkeys([name, out["aphia_name"] or name]), seen)
    if count <= 0 and out["aphia"]:
        synonyms = [s.get("scientificname") for s in
                    (_get_json(WORMS_SYN_URL.format(out["aphia"])) or [])[:25]]
        k2, c2 = _gbif_best_key([s for s in synonyms if s and SPECIES_NAME.match(s)], seen)
        if c2 > count:
            key, count = k2, c2
    if key:
        out["gbif"], out["gbif_n"] = key, max(count, 0)
    out["gbif_checklist"] = GBIF_COL_XR
    return out


def taxon_links(names, cache_file=None, online=True, workers=6, log=print):
    """
    Links for species names: {name: {name, worms, obis, gbif, worms_id, gbif_id, accepted}}.
    With online=True the WoRMS AphiaID and the GBIF species key are looked up once per name
    (kept in cache_file), giving direct links to the WoRMS page, the OBIS taxon page with its
    distribution map and the GBIF occurrence map. Without them (offline, or not in WoRMS/GBIF)
    the WoRMS and GBIF links search for the name; OBIS needs the WoRMS ID.
    """
    from urllib.parse import quote
    names = sorted({n for n in names if n})
    cache = {}
    if cache_file and os.path.exists(cache_file):
        try:
            with open(cache_file) as fh:
                cache = json.load(fh)
        except (OSError, ValueError):
            cache = {}
    todo = [n for n in names if cache.get(n, {}).get("gbif_checklist") != GBIF_COL_XR] if online else []
    if todo:
        error = None
        try:                                    # one lookup first: don't wait on every name offline
            cache[todo[0]] = _lookup_species(todo[0])
        except Exception as e:
            error = e
        if error is None and len(todo) > 1:
            with cf.ThreadPoolExecutor(max_workers=workers) as pool:
                futs = {pool.submit(_lookup_species, n): n for n in todo[1:]}
                for f in cf.as_completed(futs):
                    try:
                        cache[futs[f]] = f.result()
                    except Exception as e:
                        error = e
        missing = sum(n not in cache for n in todo)
        if missing and log:
            log(f"Taxon links: {missing} of {len(todo)} WoRMS/GBIF lookups failed ({error}); "
                f"those link to a name search instead")
        if cache_file:
            with open(cache_file, "w") as fh:
                json.dump(cache, fh, indent=1, sort_keys=True)
    links = {}
    for n in names:
        c, q = cache.get(n) or {}, quote(n)
        links[n] = {
            "name": n, "worms_id": c.get("aphia"), "gbif_id": c.get("gbif"), "gbif_n": c.get("gbif_n"),
            "worms": (f"https://www.marinespecies.org/aphia.php?p=taxdetails&id={c['aphia']}"
                      if c.get("aphia") else f"https://www.marinespecies.org/aphia.php?p=taxlist&tName={q}"),
            "obis": f"https://obis.org/taxon/{c['aphia']}" if c.get("aphia") else None,
            "gbif": (GBIF_MAP_URL.format(quote(str(c["gbif"])))
                     if c.get("gbif") else f"https://www.gbif.org/species/search?q={q}"),
            "accepted": c.get("aphia_name") if c.get("aphia_name") not in (None, n) else None}
    return links


# --------------------------------------------------------------------------
# Step 3-4: alignment and trimming
# --------------------------------------------------------------------------

def trim_to_query(aln, query_name, drop_query_gaps=False):
    q = aln[query_name]
    nongap = [i for i, c in enumerate(q) if c != "-"]
    start, end = nongap[0], nongap[-1] + 1
    cols = range(start, end)
    if drop_query_gaps:
        cols = [i for i in cols if q[i] != "-"]
    trimmed = OrderedDict()
    for name, s in aln.items():
        trimmed[name] = "".join(s[i] for i in cols)
    return trimmed, start, end


# --------------------------------------------------------------------------
# Step 5: haplotypes and minimum spanning network
# --------------------------------------------------------------------------

AMBIG = set("NRYKMSWBDHV?")


def mask_terminal_gaps(seq):
    """Leading/trailing gaps (hit shorter than the query) become missing data."""
    core = seq.strip("-")
    if not core:
        return "?" * len(seq)
    left = len(seq) - len(seq.lstrip("-"))
    right = len(seq) - len(seq.rstrip("-"))
    return "?" * left + core + "?" * right


def prepare_for_network(trimmed, gaps):
    """Return sequences used for haplotype definition and distances."""
    trimmed = OrderedDict((k, mask_terminal_gaps(v)) for k, v in trimmed.items())
    if gaps == "remove":
        L = len(next(iter(trimmed.values())))
        keep = [i for i in range(L)
                if all(s[i] not in "-?" for s in trimmed.values())]
        return OrderedDict((k, "".join(s[i] for i in keep))
                           for k, s in trimmed.items())
    return trimmed


def distance(a, b, gaps):
    d = 0
    for x, y in zip(a, b):
        if x == y or x in AMBIG or y in AMBIG:
            continue
        if gaps == "missing" and (x == "-" or y == "-"):
            continue
        d += 1
    return d


def build_haplotypes(seqs, gaps, query_name):
    """
    Collapse sequences into haplotypes. Sequences that differ only at missing
    or ambiguous positions (e.g. a shorter hit) join the same haplotype.
    Most complete sequences are processed first, so they act as representatives.
    """
    def n_missing(s):
        return sum(c in AMBIG for c in s)
    order = sorted(seqs, key=lambda k: (k != query_name, n_missing(seqs[k])))
    haps = []                                  # [rep_seq, [members]]
    for name in order:
        s = seqs[name]
        for h in haps:
            if distance(h[0], s, gaps) == 0:
                h[1].append(name)
                break
        else:
            haps.append([s, [name]])
    return [(h[0], h[1]) for h in haps]


def minimum_spanning_network(dist):
    """
    Minimum spanning network (Bandelt et al. 1999, epsilon = 0): the union
    of all minimum spanning trees.
    """
    n = len(dist)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    pairs = defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            pairs[dist[i][j]].append((i, j))
    edges = []
    for d in sorted(pairs):
        level = [(i, j) for i, j in pairs[d] if find(i) != find(j)]
        edges.extend((i, j, d) for i, j in level)
        for i, j in level:
            parent[find(i)] = find(j)
        if len({find(x) for x in range(n)}) == 1:
            break
    return edges



# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

PALETTE = ["#4e79a7", "#f28e2b", "#59a14f", "#b07aa1", "#edc948", "#76b7b2",
           "#ff9da7", "#9c755f", "#bab0ac", "#86bcb6", "#d37295", "#a0cbe8",
           "#8cd17d", "#b6992d", "#499894", "#f1ce63", "#79706e", "#d4a6c8",
           "#fabfd2", "#9d7660"]
QUERY_COLOR = "#e15759"



# --------------------------------------------------------------------------
# Per-OTU work: alignment, trimming, haplotypes, network (runs in parallel)
# --------------------------------------------------------------------------

def process_otu(job):
    """
    job: dict with otu, qname, qseq, hit_seqs [(name, seq)], outdir, mafft,
    mafft_mode, gaps, drop_query_gaps.
    Returns haplotypes as lists of member indices (-1 = query, k = k-th hit)
    and MSN edges.
    """
    d = job["outdir"]
    os.makedirs(d, exist_ok=True)
    qname = job["qname"]
    records = [(qname, job["qseq"])] + job["hit_seqs"]
    if len(records) == 1:
        return {"otu": job["otu"], "haps": [[-1]], "edges": [], "L": len(job["qseq"]),
                "span": [1, len(job["qseq"])]}

    raw = os.path.join(d, "unaligned.fasta")
    write_fasta(records, raw)
    cmd = [job["mafft"], f"--{job['mafft_mode']}", "--thread", "1", "--quiet", raw]
    if job["mafft_mode"] != "auto":
        cmd[2:2] = ["--maxiterate", "1000"]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        return {"otu": job["otu"], "error": "MAFFT failed: " + res.stderr[-500:]}
    with open(os.path.join(d, "alignment_full.fasta"), "w") as fh:
        fh.write(res.stdout)
    aln = OrderedDict((k, v.upper()) for k, v in
                      read_fasta(os.path.join(d, "alignment_full.fasta")).items())

    trimmed, start, end = trim_to_query(aln, qname, job["drop_query_gaps"])
    write_fasta(list(trimmed.items()), os.path.join(d, "alignment_trimmed.fasta"))
    L = len(next(iter(trimmed.values())))

    netseqs = prepare_for_network(trimmed, job["gaps"])
    haps = build_haplotypes(netseqs, job["gaps"], qname)
    n = len(haps)
    dist = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            dist[i][j] = dist[j][i] = distance(haps[i][0], haps[j][0], job["gaps"])
    edges = minimum_spanning_network(dist)

    index = {qname: -1}
    index.update({name: k for k, (name, _) in enumerate(job["hit_seqs"])})
    return {"otu": job["otu"],
            "haps": [[index[m] for m in members] for _, members in haps],
            "edges": [[i, j, dd] for i, j, dd in edges],
            "L": L, "span": [start + 1, end]}


# --------------------------------------------------------------------------
# HTML output
# --------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<script src="https://unpkg.com/vis-network@9.1.9/standalone/umd/vis-network.min.js"></script>
<style>
 html,body{margin:0;height:100%;font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;background:#fff;color:#222}
 #net{position:absolute;top:0;left:0;right:340px;bottom:300px}
 #empty{position:absolute;top:45%;left:0;right:340px;text-align:center;color:#888;display:none}
 #split{position:absolute;left:0;right:340px;bottom:300px;height:7px;margin-bottom:-4px;cursor:row-resize;z-index:3}
 #split:hover{background:#0001}
 #alnpane{position:absolute;left:0;right:340px;bottom:0;height:300px;border-top:1px solid #ddd;background:#fff;display:flex;flex-direction:column}
 #alnbar{display:flex;gap:8px;align-items:center;padding:5px 10px;border-bottom:1px solid #eee;background:#fafafa;font-size:12px;flex:none}
 #alnbar b{font-size:12px;text-transform:uppercase;letter-spacing:.05em;color:#666}
 #alnsel{flex:1;min-width:0;max-width:640px;font-size:12px;padding:2px}
 #alnbar button{font-size:12px}
 #aln{overflow:auto;padding:10px 14px;flex:1}
 #aln .at{font-size:16px;font-weight:600;margin-bottom:3px}
 #aln .am{color:#555;margin:2px 0}
 #aln .am a{color:#1a5fb4}
 #aln table.as{border-collapse:collapse;margin:8px 0 10px;font-size:12px}
 #aln table.as th{font-weight:normal;color:#666;text-align:left;padding:2px 22px 0 0;border-top:1px solid #ccc}
 #aln table.as td{padding:1px 22px 3px 0;border-bottom:1px solid #ccc;white-space:nowrap}
 #aln pre{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12px;line-height:1.35;margin:0}
 #aln pre b{color:#c01c28;background:#fde2e2;font-weight:600}
 #aln .none{color:#888}
 .tl a{color:#1a5fb4}
 .tl .srch{color:#888;font-size:11px}
 #side{position:absolute;top:0;right:0;width:340px;bottom:0;overflow:auto;border-left:1px solid #ddd;padding:14px;box-sizing:border-box;font-size:13px;background:#fafafa}
 h1{font-size:16px;margin:0 0 10px} h2{font-size:12px;margin:16px 0 6px;text-transform:uppercase;letter-spacing:.05em;color:#666}
 .nav{display:flex;gap:4px;align-items:center}
 select{flex:1;min-width:0;font-size:13px;padding:4px}
 .nav button{font-size:14px;width:30px}
 #pos{color:#777;font-size:12px;margin:4px 0 0}
 .lg{display:flex;align-items:center;margin:3px 0;font-style:italic}
 .lg span{display:inline-block;width:12px;height:12px;border-radius:50%;margin-right:8px;flex:none;border:1px solid #0003}
 .ctl button,.ctl label{font-size:12px;margin:2px 4px 2px 0}
 table.info{border-collapse:collapse;width:100%} table.info td{padding:2px 4px;vertical-align:top} table.info td:first-child{color:#666;white-space:nowrap}
 #detail{white-space:pre-wrap;word-break:break-word;max-height:40vh;overflow:auto;font-family:ui-monospace,Menlo,monospace;font-size:11px;background:#fff;border:1px solid #ddd;padding:8px;min-height:40px}
 input[type=search]{width:100%;box-sizing:border-box;padding:4px;margin-top:4px}
 kbd{border:1px solid #bbb;border-radius:3px;padding:0 4px;background:#fff;font-size:11px}
</style></head><body>
<div id="net"></div><div id="empty"></div>
<div id="split" title="Drag to resize the alignment pane"></div>
<div id="alnpane">
 <div id="alnbar"><b>Alignment</b><select id="alnsel" title="Top BLAST hits of this OTU"></select>
  <button id="alntoggle" title="Hide or show the alignment">Hide</button></div>
 <div id="aln"></div>
</div>
<div id="side">
 <h1>__TITLE__</h1>
 <div class="nav">
  <button id="prev" title="Previous OTU (← or ↑)">‹</button>
  <select id="otu"></select>
  <button id="next" title="Next OTU (→ or ↓)">›</button>
 </div>
 <input type="search" id="otufilter" placeholder="Filter OTU list (ID or taxon)…">
 <div id="pos"></div>
 <h2>OTU</h2>
 <table class="info" id="info"></table>
 <h2>Controls</h2>
 <div class="ctl">
 <button onclick="net.fit({animation:true})">Fit view</button>
 <button onclick="togglePhysics()" id="phys">Relax layout</button>
 <button onclick="savePNG()">Save PNG</button><br>
 <label><input type="checkbox" id="showLabels" checked onchange="render()"> labels</label>
 <label><input type="checkbox" id="showSteps" checked onchange="render()"> mutation counts</label>
 <input type="search" id="find" placeholder="Find taxon / accession in this network…" oninput="findNode(this.value)">
 </div>
 <h2>Selected haplotype</h2>
 <div id="detail">Click a node.</div>
 <h2>Legend</h2>
 <div id="legend"></div>
 <p style="color:#777;margin-top:14px"><kbd>←</kbd><kbd>→</kbd> or <kbd>↑</kbd><kbd>↓</kbd> step through OTUs.
 Scroll to zoom, drag to pan, drag nodes to rearrange. Node area ∝ number of sequences;
 edge numbers = mutational steps. Click a node to show its BLAST alignment below the network.</p>
</div>
<script>
const D = __DATA__;
const PALETTE = __PALETTE__, QCOL = "__QCOL__", GREY = "#cccccc", MAXLAB = __MAXLAB__;
const sel = document.getElementById('otu');
let visible = D.otus.map((_, i) => i), cur = 0, physicsOn = false, curNodes = [];
const esc = s => String(s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const fmt = x => (typeof x === 'number') ? x.toLocaleString() : x;

const net = new vis.Network(document.getElementById('net'), {nodes: [], edges: []}, {
  nodes: {shape:'dot', font:{size:13, face:'system-ui', strokeWidth:3, strokeColor:'#fff'}, borderWidth:1},
  edges: {color:{color:'#888', highlight:'#222'}, width:1.5, smooth:false, font:{size:11, align:'middle', background:'#fff'}},
  physics: {solver:'forceAtlas2Based', forceAtlas2Based:{gravitationalConstant:-60, springLength:60, avoidOverlap:0.6},
            stabilization:{iterations:800}},
  layout: {randomSeed: 7},
  interaction: {hover:true, tooltipDelay:100, navigationButtons:true, keyboard:false}
});
net.on('stabilizationIterationsDone', () => { net.setOptions({physics:false}); physicsOn = false;
  document.getElementById('phys').textContent = 'Relax layout'; net.fit(); });

function fillSelect() {
  sel.innerHTML = visible.map(i => { const o = D.otus[i];
    return `<option value="${i}">${esc(o.id)} — ${esc(o.assign.label)}${o.reads != null ? ' · ' + fmt(o.reads) + ' reads' : ''}</option>`; }).join('');
}

function show(i, pushHash = true) {
  cur = i; sel.value = String(i);
  const o = D.otus[i];
  const k = visible.indexOf(i);
  document.getElementById('pos').textContent = `${k + 1} of ${visible.length}` + (visible.length < D.otus.length ? ` (filtered from ${D.otus.length})` : '');
  if (pushHash) history.replaceState(null, '', '#' + encodeURIComponent(o.id));
  // info table
  const rows = [['OTU', `<b>${esc(o.id)}</b> (${o.len_full ? o.len_full + ' bp, ' + o.len + ' bp analysed' : o.len + ' bp'})`]];
  if (o.reads != null) rows.push(['Reads', fmt(o.reads)]);
  if (o.table_tax) rows.push(['Table taxonomy', esc(o.table_tax)]);
  rows.push(['Assignment', `<i>${esc(o.assign.label)}</i> (${esc(o.assign.rank || '')})`]);
  if (o.assign.n) rows.push(['Support', `${Math.round(o.assign.support * 100)}% of ${o.assign.n} best hits at ${o.assign.identity.toFixed(1)}% id`]);
  if (o.assign.taxa) rows.push(['Best-hit taxa', esc(o.assign.taxa)]);
  if (o.links) rows.push(['Best hit', `<span class="tl"><i>${esc(o.links.name)}</i>${o.links.accepted ? ' (WoRMS: <i>' + esc(o.links.accepted) + '</i>)' : ''}<br>${linkRow(o.links)}</span>`]);
  rows.push(['Network', `${o.nhits} hits, ${o.nodes.length} haplotypes, ${o.L} positions`]);
  if (o.note) rows.push(['Note', esc(o.note)]);
  document.getElementById('info').innerHTML = rows.map(r => `<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join('');
  document.getElementById('detail').textContent = 'Click a node.';
  document.getElementById('find').value = '';
  fillAln();
  render(true);
}

// ---- links to WoRMS, OBIS and GBIF for the best hit ---------------------------
function linkRow(L) {
  const a = (url, text, note) => url ? `<a href="${esc(url)}" target="_blank" rel="noopener">${text}</a>` +
    (note ? ` <span class="srch">(${note})</span>` : '') : '';
  const gbifNote = !L.gbif_id ? 'search' : L.gbif_n == null ? '' :
    L.gbif_n === 0 ? 'no records' : fmtN(L.gbif_n) + ' records';
  return [a(L.worms, 'WoRMS', L.worms_id ? '' : 'search'), a(L.obis, 'OBIS map'),
          a(L.gbif, L.gbif_id ? 'GBIF map' : 'GBIF', gbifNote)].filter(x => x).join(' · ');
}

// ---- NCBI-BLAST-style pairwise alignment -----------------------------------
let alnH = 300, alnOpen = true;
const fmtN = x => Number(x).toLocaleString('en-US');
function ncbiLink(acc) {
  const m = /^([A-Za-z]+_?\d+\.\d+)/.exec(acc);
  return m ? `<a href="https://www.ncbi.nlm.nih.gov/nuccore/${m[1]}" target="_blank" rel="noopener">${esc(acc)}</a>` : esc(acc);
}
function genbankRange(a) {           // MIDORI IDs: ACC.VER.START.END ('<' / '>' = partial)
  const m = /^[A-Za-z]+_?\d+\.\d+\.[<>]?(\d+)\.[<>]?(\d+)$/.exec(a.acc);
  if (!m) return '';
  const st = +m[1], en = +m[2], pos = p => st <= en ? st + p - 1 : st - p + 1;
  return ` &nbsp;(GenBank record: ${pos(a.ss)} to ${pos(a.se)})`;
}
function alnBlocks(a) {
  const W = 60, step = a.ss <= a.se ? 1 : -1;
  const wd = String(Math.max(a.qs, a.qe, a.ss, a.se)).length;
  const num = n => String(n).padEnd(wd);
  const out = [];
  let qi = a.qs, si = a.ss;
  for (let k = 0; k < a.q.length; k += W) {
    const qc = a.q.slice(k, k + W), sc = a.s.slice(k, k + W);
    let ql = '', sl = '', bl = '';
    for (let i = 0; i < qc.length; i++) {
      const x = qc[i], y = sc[i], ok = x === y && x !== '-';
      bl += ok ? '|' : ' ';
      ql += ok ? x : `<b>${esc(x)}</b>`;
      sl += ok ? y : `<b>${esc(y)}</b>`;
    }
    const qe = qi + qc.replace(/-/g, '').length - 1;
    const se = si + step * (sc.replace(/-/g, '').length - 1);
    out.push(`Query  ${num(qi)}  ${ql}  ${qe}`, `       ${' '.repeat(wd)}  ${bl}`,
             `Sbjct  ${num(si)}  ${sl}  ${se}`, '');
    qi = qe + 1; si = se + step;
  }
  return out.join('\n');
}
function showAln(k) {
  const o = D.otus[cur], A = o.aln || [], box = document.getElementById('aln');
  if (!A.length) {
    box.innerHTML = `<div class="none">${o.nhits ? 'No pairwise alignment stored for this OTU (only BLAST searches store one).'
                                                  : 'No BLAST hit to align.'}</div>`;
    return;
  }
  const a = A[k], pct = (x, n) => Math.round(100 * x / n);
  document.getElementById('alnsel').value = String(k);
  const cover = pct(a.qe - a.qs + 1, a.ql);
  box.innerHTML =
    `<div class="at"><i>${esc(a.label)}</i></div>` +
    (k === 0 && o.links ? `<div class="am tl">${linkRow(o.links)}</div>` : '') +
    `<div class="am">Sequence ID: ${ncbiLink(a.acc)} &nbsp; Length: ${fmtN(a.sl)} &nbsp; Number of Matches: 1</div>` +
    `<div class="am"><b>Range: ${a.ss} to ${a.se}</b>${genbankRange(a)}</div>` +
    `<table class="as"><tr><th>Score</th><th>Expect</th><th>Identities</th><th>Gaps</th><th>Strand</th><th>Query cover</th></tr>` +
    `<tr><td>${a.bits != null ? Math.round(a.bits) + ' bits' : '–'}${a.score != null ? '(' + a.score + ')' : ''}</td>` +
    `<td>${a.ev != null ? esc(a.ev) : '–'}</td><td>${a.ident}/${a.len}(${pct(a.ident, a.len)}%)</td>` +
    `<td>${a.gaps}/${a.len}(${pct(a.gaps, a.len)}%)</td><td>Plus/${a.ss <= a.se ? 'Plus' : 'Minus'}</td>` +
    `<td>${cover}% (${a.qs}–${a.qe} of ${a.ql} bp)</td></tr></table>` +
    `<pre>${alnBlocks(a)}</pre>`;
}
function fillAln() {
  const o = D.otus[cur], A = o.aln || [], sel = document.getElementById('alnsel');
  sel.innerHTML = A.map((a, k) => `<option value="${k}">${k === 0 ? 'Best hit' : 'Hit ' + (k + 1)}: ${esc(a.label)} — ${esc(a.acc)} (${a.pid.toFixed(1)}%)</option>`).join('');
  sel.disabled = !A.length;
  showAln(0);
}
function layoutAln() {
  const h = alnOpen ? alnH : 31;
  document.getElementById('alnpane').style.height = h + 'px';
  document.getElementById('net').style.bottom = h + 'px';
  document.getElementById('split').style.bottom = h + 'px';
  document.getElementById('split').style.display = alnOpen ? '' : 'none';
  document.getElementById('aln').style.display = alnOpen ? '' : 'none';
  document.getElementById('alntoggle').textContent = alnOpen ? 'Hide' : 'Show';
  net.redraw();
}
document.getElementById('alnsel').onchange = e => { showAln(Number(e.target.value)); e.target.blur(); };
document.getElementById('alntoggle').onclick = () => { alnOpen = !alnOpen; layoutAln(); };
document.getElementById('split').addEventListener('mousedown', e => {
  e.preventDefault();
  const move = ev => { alnH = Math.max(80, Math.min(window.innerHeight - 120, window.innerHeight - ev.clientY)); layoutAln(); };
  const up = () => { document.removeEventListener('mousemove', move); document.removeEventListener('mouseup', up); };
  document.addEventListener('mousemove', move); document.addEventListener('mouseup', up);
});

function render(relayout = false) {
  const o = D.otus[cur];
  const showLab = document.getElementById('showLabels').checked;
  const showSteps = document.getElementById('showSteps').checked;
  // colours by taxon, most frequent first
  const cnt = {};
  o.nodes.forEach(n => n.m.forEach(([t]) => { if (t >= 0) { const l = D.targets[t][2]; cnt[l] = (cnt[l] || 0) + 1; } }));
  const order = Object.keys(cnt).map(Number).sort((a, b) => cnt[b] - cnt[a] || D.labels[a].localeCompare(D.labels[b]));
  const col = {}; order.forEach((l, k) => col[l] = k < PALETTE.length ? PALETTE[k] : GREY);
  curNodes = o.nodes.map((n, idx) => {
    const labs = {}; let q = false, best = 0;
    n.m.forEach(([t, id]) => { if (t < 0) { q = true; return; } const l = D.targets[t][2]; labs[l] = (labs[l] || 0) + 1; best = Math.max(best, id); });
    const ls = Object.keys(labs).map(Number).sort((a, b) => labs[b] - labs[a]);
    let txt = ls.slice(0, MAXLAB).map(l => D.labels[l]).join('\n');
    if (ls.length > MAXLAB) txt += `\n+${ls.length - MAXLAB} more taxa`;
    if (q) txt = 'QUERY' + (txt ? '\n= ' + txt : '');
    const full = ls.map(l => D.labels[l] + (ls.length > 1 ? ` (${labs[l]})` : '')).join(', ');
    const det = [`Hap_${idx + 1}  (n = ${n.m.length})`];
    if (q) det.push('Contains the query (OTU) sequence');
    if (ls.length) { det.push(`Identity to query: ${best.toFixed(1)}%`); det.push('Taxa: ' + full); }
    det.push('Members:');
    n.m.forEach(([t, id]) => det.push(t < 0 ? `  ${o.id} (query)` : `  ${D.targets[t][0]}  ${id.toFixed(1)}%  ${D.targets[t][1]}`));
    return {id: idx, label: showLab ? txt : ' ', size: 8 + 6 * Math.sqrt(n.m.length),
      color: q ? {background: QCOL, border: '#7a1416'} : {background: col[ls[0]], border: '#333'},
      font: q ? {size: 16, bold: true} : {size: 13},
      title: `Hap_${idx + 1} (n=${n.m.length})${q ? ' — QUERY' : ''}: ${full}`.slice(0, 300),
      detail: det.join('\n')};
  });
  const edges = o.edges.map(([a, b, d], k) => ({id: k, from: a, to: b, label: showSteps ? String(d) : ' ',
    length: 40 + 25 * Math.min(d, 30), title: `${d} mutational step(s)`}));
  if (relayout) {
    net.setData({nodes: new vis.DataSet(curNodes), edges: new vis.DataSet(edges)});
    net.setOptions({physics: true}); physicsOn = true;
    document.getElementById('phys').textContent = 'Freeze layout';
    document.getElementById('empty').style.display = 'none';
  } else {
    net.body.data.nodes.update(curNodes.map(n => ({id: n.id, label: n.label})));
    net.body.data.edges.update(edges.map(e => ({id: e.id, label: e.label})));
  }
  // legend
  const lg = [`<div class="lg"><span style="background:${QCOL}"></span><b style="font-style:normal">Query</b></div>`];
  order.slice(0, PALETTE.length).forEach(l => lg.push(`<div class="lg"><span style="background:${col[l]}"></span>${esc(D.labels[l])} (${cnt[l]})</div>`));
  if (order.length > PALETTE.length) lg.push(`<div class="lg"><span style="background:${GREY}"></span><span style="all:unset;font-style:normal">${order.length - PALETTE.length} other taxa</span></div>`);
  document.getElementById('legend').innerHTML = lg.join('');
}

function step(delta) {
  if (!visible.length) return;
  let k = visible.indexOf(cur);
  k = k < 0 ? 0 : Math.min(visible.length - 1, Math.max(0, k + delta));
  show(visible[k]);
}
document.getElementById('prev').onclick = () => step(-1);
document.getElementById('next').onclick = () => step(1);
sel.onchange = () => { show(Number(sel.value)); sel.blur(); };
document.addEventListener('keydown', e => {
  const tag = (e.target.tagName || '').toLowerCase();
  if (tag === 'input' || tag === 'textarea' || tag === 'select' || e.altKey || e.ctrlKey || e.metaKey) return;
  if (e.key === 'ArrowRight' || e.key === 'ArrowDown') { step(1); e.preventDefault(); }
  else if (e.key === 'ArrowLeft' || e.key === 'ArrowUp') { step(-1); e.preventDefault(); }
  else if (e.key === 'Home') { step(-1e9); e.preventDefault(); }
  else if (e.key === 'End') { step(1e9); e.preventDefault(); }
});
document.getElementById('otufilter').addEventListener('input', e => {
  const q = e.target.value.trim().toLowerCase();
  visible = D.otus.map((o, i) => i).filter(i => { const o = D.otus[i];
    return !q || (o.id + ' ' + o.assign.label + ' ' + (o.table_tax || '') + ' ' + (o.assign.taxa || '')).toLowerCase().includes(q); });
  fillSelect();
  if (visible.length && !visible.includes(cur)) show(visible[0]); else if (visible.length) show(cur, false);
});
function togglePhysics() { physicsOn = !physicsOn; net.setOptions({physics: physicsOn});
  document.getElementById('phys').textContent = physicsOn ? 'Freeze layout' : 'Relax layout'; }
function savePNG() { const c = document.querySelector('#net canvas'), a = document.createElement('a');
  a.download = `haplotype_network_${D.otus[cur].id}.png`; a.href = c.toDataURL('image/png'); a.click(); }
net.on('click', p => {
  if (!p.nodes.length) return;
  document.getElementById('detail').textContent = curNodes[p.nodes[0]].detail;
  const o = D.otus[cur], A = o.aln || [];
  const members = new Set(o.nodes[p.nodes[0]].m.map(([t]) => t));
  const k = A.findIndex(a => members.has(a.t));
  if (k >= 0) { showAln(k); if (!alnOpen) { alnOpen = true; layoutAln(); } }
});
function findNode(q) { q = q.trim().toLowerCase(); if (!q) { net.unselectAll(); return; }
  const hits = curNodes.filter(n => n.detail.toLowerCase().includes(q)).map(n => n.id);
  net.selectNodes(hits); if (hits.length) net.fit({nodes: hits, animation: true}); }

fillSelect();
const start = D.otus.findIndex(o => '#' + encodeURIComponent(o.id) === location.hash);
show(start >= 0 ? start : 0, false);
</script></body></html>"""


def write_html(path, title, data, max_label_taxa):
    page = (HTML_TEMPLATE
            .replace("__TITLE__", html.escape(title))
            .replace("__PALETTE__", json.dumps(PALETTE))
            .replace("__QCOL__", QUERY_COLOR)
            .replace("__MAXLAB__", str(int(max_label_taxa)))
            .replace("__DATA__", json.dumps(data, separators=(",", ":"))
                     .replace("</", "<\\/")))
    with open(path, "w") as fh:
        fh.write(page)


# --------------------------------------------------------------------------
# Input
# --------------------------------------------------------------------------

def read_input(path_or_seq, id_col=None, seq_col=None):
    """
    Return (rows, fieldnames, delimiter, id_col, seq_col) where rows is a list
    of dicts. For FASTA / raw-sequence input fieldnames is None.
    """
    if not os.path.isfile(path_or_seq):
        if not re.fullmatch(r"[ACGTUNRYKMSWBDHVacgtunrykmswbdhv\s\-]+", path_or_seq):
            sys.exit(f"ERROR: input file '{path_or_seq}' not found (working directory: "
                     f"{os.getcwd()}). Give the correct path, or a raw DNA sequence.")
        seq = re.sub(r"[^A-Za-z]", "", path_or_seq)
        return [{"OTU": "query", "sequence": seq}], None, None, "OTU", "sequence"
    with open(path_or_seq) as fh:
        first = fh.readline()
    if first.startswith(">"):
        rows = [{"OTU": seq_id(h), "sequence": s} for h, s in read_fasta(path_or_seq).items()]
        return rows, None, None, "OTU", "sequence"

    delim = "\t" if first.count("\t") > first.count(",") else ","
    with open(path_or_seq, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        fields = reader.fieldnames
        rows = list(reader)
    low = {f.lower().strip(): f for f in fields}
    if not id_col:
        id_col = next((low[k] for k in ("otu", "#otu id", "otu_id", "otuid", "id", "asv",
                                         "asv_id", "zotu") if k in low), fields[0])
    if not seq_col:
        seq_col = next((low[k] for k in ("sequence", "seq", "sequences", "dna", "repseq")
                        if k in low), None)
    if not seq_col or seq_col not in fields:
        sys.exit("ERROR: no sequence column found in the table; name it with --seq-col")
    if id_col not in fields:
        sys.exit(f"ERROR: ID column '{id_col}' not in the table")
    return rows, fields, delim, id_col, seq_col


def table_taxonomy(row):
    """Lowest non-empty rank from domain..species columns of the input table."""
    low = {k.lower(): k for k in row}
    lowest = None
    for rank in RANK_ORDER:
        k = low.get(rank)
        if k and row[k].strip():
            lowest = (row[k].strip(), rank)
    if not lowest:
        rk = low.get("assigned_rank")
        return "unassigned" if rk and row[rk].strip() == "unassigned" else None
    txt = f"{lowest[0]} ({lowest[1]})"
    for c in ("identity", "pident"):
        if c in low and row[low[c]].strip():
            txt += f", {row[low[c]].strip()}% id"
            break
    return txt


def read_count(row, fields):
    low = {k.lower(): k for k in row}
    for c in ("total_reads", "total", "reads", "abundance", "size", "count"):
        if c in low:
            try:
                return int(float(row[low[c]]))
            except ValueError:
                pass
    return None


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", "-q", "--query", dest="input", required=True,
                    help="OTU table (CSV/TSV), FASTA file, or one raw DNA sequence")
    ap.add_argument("-d", "--db", required=True, help="reference database FASTA")
    ap.add_argument("-o", "--outdir", default="haplonet_out", help="output directory")
    ap.add_argument("--id-col", help="OTU ID column of the table (auto-detected)")
    ap.add_argument("--seq-col", help="sequence column of the table (auto-detected)")
    ap.add_argument("--otus", help="only process these OTU IDs (comma-separated)")
    ap.add_argument("--min-reads", type=int, default=0,
                    help="skip OTUs with fewer total reads (needs a read-count column)")
    ap.add_argument("--taxonomy", help="optional TSV: reference ID <tab> taxonomy string")
    ap.add_argument("--search", choices=["blast", "vsearch"], default="blast",
                    help="similarity search tool (default: blast)")
    ap.add_argument("--blastdb",
                    help="existing BLAST database prefix built from --db (made with "
                         "makeblastdb, without -parse_seqids); built automatically if omitted")
    ap.add_argument("--blast-task", default="megablast",
                    choices=["megablast", "dc-megablast", "blastn"],
                    help="blastn task (default megablast, as on the NCBI website; "
                         "use dc-megablast or blastn to find more divergent hits)")
    ap.add_argument("--no-trim-query", action="store_true",
                    help="(blast) do not trim OTUs to the part covered by their best hit")
    ap.add_argument("--trim-min", type=int, default=10,
                    help="(blast) trim only if at least this many bases are uncovered (10)")
    ap.add_argument("--evalue", type=float, default=1e-5, help="BLAST e-value cutoff (1e-5)")
    ap.add_argument("--min-qcov", type=float, default=80,
                    help="minimum query coverage (%%) of a hit (default 80)")
    ap.add_argument("--top", type=int, default=100, help="number of best hits to keep (100)")
    ap.add_argument("--min-id", type=float, default=0.70,
                    help="minimum identity of hits, 0-1 (0.70)")
    ap.add_argument("--max-ambig", type=float, default=0.02,
                    help="discard hits whose aligned reference region has more than this "
                         "fraction of N/IUPAC ambiguity codes (default 0.02). VSEARCH "
                         "counts N as a match, so runs of N give false 100%% hits")
    ap.add_argument("--oversample", type=int, default=5,
                    help="request top x this many hits from the search so that enough "
                         "remain after discarding ambiguous ones (5)")
    ap.add_argument("--maxrejects", type=int, default=0,
                    help="(vsearch only) --maxrejects; 0 = exhaustive search so the best hits "
                         "are really the best (slower on huge databases)")
    ap.add_argument("--strand", choices=["plus", "both"], default="both",
                    help="(vsearch only; BLAST always searches both strands)")
    ap.add_argument("--tie-delta", type=float, default=0.0,
                    help="hits within this many %%-points of the best are pooled "
                         "for the consensus assignment (0 = exact ties only)")
    ap.add_argument("--lca-cutoff", type=float, default=0.6,
                    help="fraction of the best hits that must agree at a rank "
                         "(1.0 = strict LCA; default 0.6)")
    ap.add_argument("--full-hits", action="store_true",
                    help="align full reference sequences instead of only the region "
                         "VSEARCH aligned to the OTU")
    ap.add_argument("--max-label-taxa", type=int, default=2,
                    help="taxa named in a node label before '+N more' (2)")
    ap.add_argument("--no-taxon-links", action="store_true",
                    help="don't look up WoRMS/GBIF IDs online for the best-hit links "
                         "(the links then search for the species name)")
    ap.add_argument("--max-alignments", type=int, default=10,
                    help="top hits whose pairwise BLAST alignment is shown in the page (10)")
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--mafft-mode", default="auto",
                    choices=["auto", "localpair", "globalpair", "genafpair"])
    ap.add_argument("--drop-query-gaps", action="store_true",
                    help="also remove columns where the OTU has a gap")
    ap.add_argument("--gaps", choices=["fifth", "missing", "remove"], default="fifth",
                    help="gap treatment in the network (default: fifth state)")
    ap.add_argument("--title", default="Haplotype networks", help="page title")
    ap.add_argument("--vsearch", default="vsearch", help="path to vsearch")
    ap.add_argument("--blastn", default="blastn", help="path to blastn")
    ap.add_argument("--makeblastdb", default="makeblastdb", help="path to makeblastdb")
    ap.add_argument("--mafft", default="mafft", help="path to mafft")
    args = ap.parse_args()

    if not os.path.isfile(args.db):
        sys.exit(f"ERROR: reference FASTA '{args.db}' not found")
    if args.search == "blast" and args.blastdb:
        prefix = args.blastdb
        for ext in (".nsq", ".nal", ".ndb", ".nin"):
            if prefix.endswith(ext):          # accept 'MIDORI2_lrRNA.nsq' too
                prefix = prefix[:-len(ext)]
        if not any(os.path.exists(prefix + e) for e in (".nal", ".nin", ".00.nin", ".ndb")):
            d = os.path.dirname(prefix) or "."
            found = sorted({f.rsplit(".", 1)[0] for f in os.listdir(d)
                            if f.endswith((".nin", ".nal"))}) if os.path.isdir(d) else []
            hint = (" BLAST databases in that folder: " +
                    ", ".join(os.path.join(d, f) if d != "." else f for f in found)
                    if found else "")
            sys.exit(f"ERROR: no BLAST database '{prefix}' found.{hint}\n"
                     f"Give --blastdb the database prefix (the file name without the "
                     f".nsq/.nin/... extension), or leave --blastdb out to build one.")
        args.blastdb = prefix
    if args.search == "blast":
        blastn = check_tool("blast", args.blastn)
        makeblastdb = check_tool("blast", args.makeblastdb) if not args.blastdb else None
    else:
        vsearch = check_tool("vsearch", args.vsearch)
    mafft = check_tool("mafft", args.mafft)
    os.makedirs(args.outdir, exist_ok=True)
    log = os.path.join(args.outdir, "tools.log")
    open(log, "w").close()

    # ---- Input -----------------------------------------------------------
    rows, fields, delim, id_col, seq_col = read_input(args.input, args.id_col, args.seq_col)
    if args.otus:
        wanted = {x.strip() for x in args.otus.split(",")}
        rows = [r for r in rows if r[id_col] in wanted]
    if args.min_reads:
        rows = [r for r in rows if (read_count(r, fields) or 0) >= args.min_reads]
    otus, used = [], set()
    for r in rows:
        seq = re.sub(r"[^A-Za-z]", "", r[seq_col]).upper().replace("U", "T")
        if not seq:
            continue
        safe = re.sub(r"[^\w.\-]", "_", r[id_col]) or "otu"
        base, k = safe, 2
        while safe in used:
            safe, k = f"{base}_{k}", k + 1
        used.add(safe)
        otus.append({"id": r[id_col], "safe": safe, "seq": seq, "row": r})
    if not otus:
        sys.exit("ERROR: no sequences to process")
    print(f"[1/5] {len(otus)} sequences to process", file=sys.stderr)

    query_fa = os.path.join(args.outdir, "queries.fasta")
    write_fasta([(o["safe"], o["seq"]) for o in otus], query_fa)

    # ---- Reference & taxonomy -------------------------------------------
    print("[2/5] Reading reference database", file=sys.stderr)
    ref = read_fasta(args.db)
    ref_by_id = {}
    for h, s in ref.items():
        ref_by_id[h] = (h, s)
        ref_by_id.setdefault(seq_id(h), (h, s))
        ref_by_id.setdefault(h.replace("###", " ", 1).split()[0], (h, s))   # BLAST sseqid
    tax_tsv = load_taxonomy_tsv(args.taxonomy) if args.taxonomy else {}
    tax_cache = {}

    def tax_of_target(target):
        if target not in tax_cache:
            header = ref_by_id[target][0]
            tstr = tax_tsv.get(seq_id(header)) or taxonomy_from_header(header)
            lineage = parse_taxonomy_string(tstr)
            label, rank = lowest_label(lineage)
            tax_cache[target] = {"lineage": lineage, "label": label, "rank": rank}
        return tax_cache[target]

    # ---- Similarity search (all OTUs in one run) ------------------------
    n_request = args.top * max(1, args.oversample)
    if args.search == "blast":
        blastdb = args.blastdb
        if not blastdb:
            print("[3/5] Building BLAST database", file=sys.stderr)
            blastdb = make_blast_db(makeblastdb, args.db, args.outdir, log)
        print(f"[3/5] BLAST search ({args.blast_task}, all sequences)", file=sys.stderr)
        by_query = blast_hits(blastn, query_fa, blastdb, args.outdir, n_request,
                              args.min_id, args.threads, args.blast_task,
                              args.evalue, args.min_qcov, log)
        n_trimmed = 0
        for o in otus:
            hits, lo, hi = fit_hits_to_query(by_query.get(o["safe"], []), len(o["seq"]),
                                             args.min_qcov, not args.no_trim_query,
                                             args.trim_min)
            by_query[o["safe"]] = hits
            if (lo, hi) != (1, len(o["seq"])):
                o["trim"] = (lo, hi, len(o["seq"]))
                o["seq"] = o["seq"][lo - 1:hi]
                n_trimmed += 1
        if n_trimmed:
            print(f"      {n_trimmed} sequences trimmed to the part covered by their best "
                  f"BLAST hit (see hn_query_region)", file=sys.stderr)
        missing = {h["target"] for hs in by_query.values() for h in hs} - set(ref_by_id)
        if missing:
            sys.exit(f"ERROR: {len(missing)} BLAST hit IDs (e.g. '{next(iter(missing))}') "
                     f"are not in {args.db}. Was --blastdb built from this FASTA?")
    else:
        print("[3/5] VSEARCH global search (all sequences)", file=sys.stderr)
        by_query = vsearch_hits(vsearch, query_fa, args.db, args.outdir, n_request,
                                args.min_id, args.threads, args.strand,
                                args.maxrejects, log)
        by_query = {q: [h for h in hs if float(h["qcov"]) >= args.min_qcov]
                    for q, hs in by_query.items()}

    # ---- Drop hits on N-runs / ambiguous regions, keep the best `top` ----
    # VSEARCH counts N (and other IUPAC codes) as matches, so a reference region
    # that is a run of Ns "matches" any query at 100% identity. BLAST does not,
    # but the filter is kept as a safety net for both.
    def hit_segment(h):
        s = ref_by_id[h["target"]][1]
        if not args.full_hits and h.get("tilo") and h.get("tihi"):
            s = s[int(h["tilo"]) - 1:int(h["tihi"])]
        return s.upper().replace("U", "T")

    n_dropped = {}
    fields_out = ["query", "target", "id", "alnlen", "mism", "opens", "qcov", "tcov",
                  "ql", "tl", "qstrand", "tilo", "tihi", "ambig_frac"]
    hits_file = "blast_hits.tsv" if args.search == "blast" else "vsearch_hits.tsv"
    if args.search == "blast":
        fields_out += ["evalue", "bitscore"]
    with open(os.path.join(args.outdir, hits_file), "w") as fh:
        fh.write("\t".join(fields_out) + "\n")
        for q, hits in by_query.items():
            kept, dropped = [], 0
            for h in hits:
                seg = hit_segment(h)
                h["seg"] = seg
                h["ambig_frac"] = round(sum(c not in "ACGT" for c in seg) / max(1, len(seg)), 4)
                if h["ambig_frac"] > args.max_ambig:
                    dropped += 1
                else:
                    kept.append(h)
            by_query[q] = kept[:args.top]
            n_dropped[q] = dropped
            for h in by_query[q]:
                fh.write("\t".join(str(h[f]) for f in fields_out) + "\n")
    tot_drop = sum(n_dropped.values())
    if tot_drop:
        print(f"      Discarded {tot_drop} hits whose reference region is more than "
              f"{args.max_ambig:.0%} N/ambiguous ({sum(1 for v in n_dropped.values() if v)}"
              f" sequences affected)", file=sys.stderr)

    jobs = []
    for o in otus:
        hits = by_query.get(o["safe"], [])
        o["hits"] = hits
        tax_of = {h["target"]: tax_of_target(h["target"]) for h in hits}
        o["assign"] = consensus_assignment(hits, tax_of, args.tie_delta, args.lca_cutoff)
        hit_seqs = []
        for k, h in enumerate(hits):
            s = h["seg"]
            if h.get("qstrand") == "-":
                s = revcomp(s)
            hit_seqs.append((f"H{k + 1:03d}_" + re.sub(r"[^\w.\-]", "_",
                             seq_id(h["target"]))[:60], s.upper().replace("U", "T")))
        jobs.append({"otu": o["safe"], "qname": "QUERY_" + o["safe"], "qseq": o["seq"],
                     "hit_seqs": hit_seqs,
                     "outdir": os.path.join(args.outdir, "otus", o["safe"]),
                     "mafft": mafft, "mafft_mode": args.mafft_mode, "gaps": args.gaps,
                     "drop_query_gaps": args.drop_query_gaps})

    # ---- Alignment + networks (parallel) --------------------------------
    print(f"[4/5] Aligning and building networks ({args.threads} workers)", file=sys.stderr)
    results = {}
    with cf.ProcessPoolExecutor(max_workers=args.threads) as ex:
        futs = {ex.submit(process_otu, j): j["otu"] for j in jobs}
        for k, f in enumerate(cf.as_completed(futs), 1):
            try:
                results[futs[f]] = f.result()
            except Exception as e:                       # keep going
                results[futs[f]] = {"otu": futs[f], "error": repr(e)}
            if k % 25 == 0 or k == len(futs):
                print(f"      {k}/{len(futs)} done", file=sys.stderr)

    # ---- Collect, write tables and HTML ---------------------------------
    print("[5/5] Writing outputs", file=sys.stderr)
    targets, t_index, labels, l_index = [], {}, [], {}

    def tidx(target):
        if target not in t_index:
            t = tax_of_target(target)
            lab = t["label"]
            if lab not in l_index:
                l_index[lab] = len(labels)
                labels.append(lab)
            t_index[target] = len(targets)
            targets.append([seq_id(target), "; ".join(n for r, n in t["lineage"] if r is None or r in RANK_ORDER),
                            l_index[lab]])
        return t_index[target]

    # links to WoRMS, OBIS and GBIF for best hits that are eukaryote species
    link_name = {}
    for o in otus:
        bh = best_hit(o["hits"])
        link_name[o["safe"]] = link_species(tax_of_target(bh["target"])) if bh else None
    links = taxon_links(link_name.values(), os.path.join(args.outdir, "taxon_links.json"),
                        online=not args.no_taxon_links)

    page_otus, summary = [], []
    for o in otus:
        res = results[o["safe"]]
        hits, a = o["hits"], o["assign"]
        note = None
        if not hits:
            note = f"No usable hits >= {args.min_id:.0%} identity"
        if o.get("trim"):
            lo, hi, L0 = o["trim"]
            cut = ", ".join(x for x in (f"1–{lo - 1}" if lo > 1 else "",
                                        f"{hi + 1}–{L0}" if hi < L0 else "") if x)
            note = ((note + "; ") if note else "") + (
                f"Bases {cut} of {L0} not matched by the best hit were removed; "
                f"analysed region {lo}–{hi} (possible primer artefact or chimera)")
        if n_dropped.get(o["safe"]):
            note = ((note + "; ") if note else "") + (
                f"{n_dropped[o['safe']]} hits on N-runs/ambiguous reference regions discarded")
        if "error" in res:
            note = res["error"]
            res = {"haps": [[-1] + list(range(len(hits)))], "edges": [], "L": len(o["seq"])}
        nodes = []
        for members in res["haps"]:
            m = [[-1, 100.0] if x < 0 else [tidx(hits[x]["target"]), hits[x]["id"]]
                 for x in members]
            nodes.append({"m": m})
        taxa_txt = ", ".join(f"{k} ({v})" for k, v in a["top_labels"].most_common(6)) if a else ""
        if a and len(a["top_labels"]) > 6:
            taxa_txt += ", …"
        aln = alignment_records(hits, lambda t: tax_of_target(t)["label"], args.max_alignments)
        for rec in aln:
            rec["t"] = tidx(rec.pop("target"))
        page_otus.append({
            "id": o["id"], "len": len(o["seq"]), "aln": aln,
            "len_full": o["trim"][2] if o.get("trim") else None,
            "reads": read_count(o["row"], fields) if fields else None,
            "table_tax": table_taxonomy(o["row"]) if fields else None,
            "assign": ({"label": a["label"], "rank": a["rank"], "support": a["support"],
                        "identity": a["identity"], "n": a["n_top_hits"], "taxa": taxa_txt}
                       if a else {"label": "no hits", "rank": ""}),
            "nhits": len(hits), "L": res["L"], "nodes": nodes, "edges": res["edges"],
            "links": links.get(link_name[o["safe"]]), "note": note})
        summary.append({
            "hn_assignment": a["label"] if a else "no hits",
            "hn_rank": (a["rank"] or "") if a else "",
            "hn_support": f"{a['support']:.2f}" if a else "",
            "hn_best_identity": f"{a['identity']:.2f}" if a else "",
            "hn_n_best_hits": a["n_top_hits"] if a else 0,
            "hn_taxa_best_hits": "; ".join(f"{k} ({v})" for k, v in
                                           a["top_labels"].most_common()) if a else "",
            "hn_n_hits": len(hits), "hn_n_haplotypes": len(nodes),
            "hn_n_ambiguous_hits_dropped": n_dropped.get(o["safe"], 0),
            "hn_query_region": (f"{o['trim'][0]}-{o['trim'][1]}" if o.get("trim")
                                else f"1-{len(o['seq'])}"),
            "hn_trimmed_bp": (o["trim"][2] - len(o["seq"])) if o.get("trim") else 0})

    # assignment table (input table + new columns, or a plain table)
    out_cols = list(summary[0].keys())
    if fields:
        out_path = os.path.join(args.outdir, "otu_table_haplonet" +
                                (".tsv" if delim == "\t" else ".csv"))
        with open(out_path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter=delim)
            w.writerow(fields + out_cols)
            for o, s in zip(otus, summary):
                w.writerow([o["row"][f] for f in fields] + [s[c] for c in out_cols])
    out_tsv = os.path.join(args.outdir, "taxonomic_assignments.tsv")
    with open(out_tsv, "w") as fh:
        fh.write("otu\t" + "\t".join(out_cols) + "\n")
        for o, s in zip(otus, summary):
            fh.write(o["id"] + "\t" + "\t".join(str(s[c]) for c in out_cols) + "\n")

    html_path = os.path.join(args.outdir, "haplotype_networks.html")
    write_html(html_path, args.title,
               {"otus": page_otus, "targets": targets, "labels": labels},
               args.max_label_taxa)

    n_err = sum(1 for r in results.values() if "error" in r)
    print(f"\nDone: {len(otus)} sequences{f', {n_err} with errors' if n_err else ''}.\n"
          f"Outputs in {args.outdir}/:\n"
          f"  haplotype_networks.html     interactive networks (dropdown / arrow keys)\n"
          f"  taxonomic_assignments.tsv   assignment per OTU"
          + ("\n  otu_table_haplonet.*        input table + assignment columns" if fields else "")
          + f"\n  {hits_file:<27} best hits of every OTU\n"
          f"  otus/<OTU>/                 alignments per OTU", file=sys.stderr)


if __name__ == "__main__":
    main()
