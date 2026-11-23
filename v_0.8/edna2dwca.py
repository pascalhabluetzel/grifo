#!/usr/bin/env python3
"""
edna2dwca.py  -  Metabarcoding OTU/ASV table + FAIRe-style sample metadata
                 -> zipped Darwin Core Archive (Occurrence core + DNA-derived data extension)

Follows the data model recommended for eDNA in GBIF/OBIS ("Publishing DNA-derived data
through biodiversity data platforms") and referenced by eDNAqua-Plan D3.2 (FAIRe checklist,
DNA-derived extension <-> MIxS alignment):

  * one occurrence per (sample x OTU) with reads > 0, basisOfRecord = MaterialSample
  * organismQuantity = reads, sampleSizeValue = total reads in the sample (both "DNA sequence reads")
  * taxonomy parsed from the classifier lineage; NCBI taxid kept as taxonID
  * MIxS / FAIRe assay fields go into the DNA-derived data extension
  * a QC report lists everything that is missing, inferred or suspicious

Standard library only (Python >= 3.8).

GRIFO copy
----------
Copied from MarGenoDaT (https://github.com/pascalhabluetzel/MarGenoDaT, tools/edna2dwca.py,
commit 4baab80, MIT licence) and used by section 7 of grifo.ipynb. Changes to the original,
all backwards compatible:
  * sniff_read() also reads semicolon-separated files and skips the '#' guidance rows that the
    metadata template builder can add under the header.
  * The FAIRe column samp_collect_method is used for samp_collec_method when that is empty.
  * --identity-label sets how the identity column is described in identificationRemarks
    (default unchanged: "mean identity to references"); GRIFO passes the best-hit identity.
  * All terms of the GBIF DNA-derived data extension that apply to metabarcoding are written
    when a value is available (sample metadata, config "defaults", config "sample_values" or
    an OTU-table column of that name), not only a fixed subset.
  * FAIRe checklist terms that are neither Darwin Core nor DNA-derived extension terms go into
    dynamicProperties (JSON) under their FAIRe names; the OTU ID is written there as seq_id
    (was otu_id).
  * The OTU table may carry per-OTU term columns (e.g. percent_match, scientificNameID,
    ampliconSize, verbatimIdentification); they are not mistaken for samples.
  * Occurrence core: scientificNameID and associatedSequences added.
  * Config "sample_values": {sample: {term: value}} fills per-sample values computed by the
    pipeline (read counts, ...) where the metadata has none.

Usage
-----
  python edna2dwca.py --otu otus.tsv --metadata metadata.csv -o my_dataset
  python edna2dwca.py --otu otus.tsv --metadata metadata.csv --config dataset.json \
                      --fasta otus.fasta -o my_dataset
  python edna2dwca.py --write-config-template dataset.json     # starter config to fill in

Output: my_dataset.zip (occurrence.txt, dna_derived_data.txt, meta.xml, eml.xml)
        my_dataset_QC_report.md
"""

import argparse
import csv
import datetime as dt
import io
import json
import re
import sys
import uuid
import zipfile
from collections import Counter, OrderedDict
from xml.sax.saxutils import escape

# --------------------------------------------------------------------------------------
# Term definitions
# --------------------------------------------------------------------------------------
DWC = "http://rs.tdwg.org/dwc/terms/"
DC = "http://purl.org/dc/terms/"

OCCURRENCE_TERMS = OrderedDict([
    ("occurrenceID", DWC), ("basisOfRecord", DWC), ("occurrenceStatus", DWC),
    ("eventID", DWC), ("parentEventID", DWC), ("eventDate", DWC), ("verbatimEventDate", DWC),
    ("materialSampleID", DWC), ("recordedBy", DWC), ("samplingProtocol", DWC),
    ("locality", DWC), ("waterBody", DWC), ("countryCode", DWC),
    ("decimalLatitude", DWC), ("decimalLongitude", DWC), ("geodeticDatum", DWC),
    ("coordinateUncertaintyInMeters", DWC),
    ("minimumDepthInMeters", DWC), ("maximumDepthInMeters", DWC),
    ("associatedSequences", DWC),
    ("scientificName", DWC), ("scientificNameID", DWC), ("taxonRank", DWC), ("taxonID", DWC),
    ("nameAccordingTo", DWC),
    ("identificationQualifier", DWC),
    ("kingdom", DWC), ("phylum", DWC), ("class", DWC), ("order", DWC),
    ("superfamily", DWC), ("family", DWC), ("subfamily", DWC), ("tribe", DWC), ("subtribe", DWC),
    ("genus", DWC), ("genericName", DWC), ("subgenus", DWC), ("infragenericEpithet", DWC),
    ("specificEpithet", DWC), ("infraspecificEpithet", DWC), ("higherClassification", DWC),
    ("verbatimIdentification", DWC), ("identificationRemarks", DWC),
    ("identificationReferences", DWC),
    ("organismQuantity", DWC), ("organismQuantityType", DWC),
    ("sampleSizeValue", DWC), ("sampleSizeUnit", DWC),
    ("datasetName", DWC), ("institutionCode", DWC),
    ("dynamicProperties", DWC),
])
OCC_URI = {"rowType": DWC + "Occurrence"}

MIXS = "https://w3id.org/mixs/"
GBIFT = "http://rs.gbif.org/terms/"
MIQE = "http://rs.gbif.org/terms/miqe/"
GGBN = "http://data.ggbn.org/schemas/ggbn/terms/"
# Terms of the GBIF DNA-derived data extension (http://rs.gbif.org/terms/1.0/DNADerivedData,
# issued 2024-07-11) that apply to amplicon metabarcoding, in the order of the extension.
# Genome-assembly, single-cell, virus and qPCR/ddPCR-only terms are left out.
DNA_TERMS = OrderedDict([
    ("samp_name", MIXS + "0001107"),
    ("project_name", MIXS + "0000092"),
    ("experimental_factor", MIXS + "0000008"),
    ("samp_taxon_id", MIXS + "0001320"),
    ("neg_cont_type", MIXS + "0001321"),
    ("pos_cont_type", MIXS + "0001322"),
    ("env_broad_scale", MIXS + "0000012"),
    ("env_local_scale", MIXS + "0000013"),
    ("env_medium", MIXS + "0000014"),
    ("source_mat_id", MIXS + "0000026"),
    ("samp_collec_device", MIXS + "0000002"),
    ("samp_collec_method", MIXS + "0001225"),
    ("samp_mat_process", MIXS + "0000016"),
    ("size_frac", MIXS + "0000017"),
    ("samp_size", MIXS + "0000001"),
    ("samp_vol_we_dna_ext", MIXS + "0000111"),
    ("nucl_acid_ext", MIXS + "0000037"),
    ("nucl_acid_amp", MIXS + "0000038"),
    ("lib_size", MIXS + "0000039"),
    ("lib_reads_seqd", MIXS + "0000040"),
    ("lib_layout", MIXS + "0000041"),
    ("lib_vector", MIXS + "0000042"),
    ("lib_screen", MIXS + "0000043"),
    ("target_gene", MIXS + "0000044"),
    ("target_subfragment", MIXS + "0000045"),
    ("pcr_primers", MIXS + "0000046"),
    ("mid", MIXS + "0000047"),
    ("adapters", MIXS + "0000048"),
    ("pcr_cond", MIXS + "0000049"),
    ("seq_meth", MIXS + "0000050"),
    ("seq_quality_check", MIXS + "0000051"),
    ("chimera_check", MIXS + "0000052"),
    ("tax_ident", MIXS + "0000053"),
    ("ref_db", MIXS + "0000062"),
    ("sim_search_meth", MIXS + "0000063"),
    ("tax_class", MIXS + "0000064"),
    ("contam_screen_input", MIXS + "0000005"),
    ("contam_screen_param", MIXS + "0000073"),
    ("decontam_software", MIXS + "0000074"),
    ("otu_class_appr", MIXS + "0000085"),
    ("otu_seq_comp_appr", MIXS + "0000086"),
    ("otu_db", MIXS + "0000087"),
    ("url", MIXS + "0000091"),
    ("sop", MIXS + "0000090"),
    ("pcr_primer_forward", GBIFT + "pcr_primer_forward"),
    ("pcr_primer_reverse", GBIFT + "pcr_primer_reverse"),
    ("pcr_primer_name_forward", GBIFT + "pcr_primer_name_forward"),
    ("pcr_primer_name_reverse", GBIFT + "pcr_primer_name_reverse"),
    ("pcr_primer_reference", GBIFT + "pcr_primer_reference"),
    ("DNA_sequence", GBIFT + "dna_sequence"),
    ("concentration", GGBN + "concentration"),
    ("concentrationUnit", GGBN + "concentrationUnit"),
    ("methodDeterminationConcentrationAndRatios", GGBN + "methodDeterminationConcentrationAndRatios"),
    ("ratioOfAbsorbance260_230", GGBN + "ratioOfAbsorbance260_230"),
    ("ratioOfAbsorbance260_280", GGBN + "ratioOfAbsorbance260_280"),
    ("annealingTemp", MIQE + "annealingTemp"),
    ("annealingTempUnit", MIQE + "annealingTempUnit"),
    ("ampliconSize", MIQE + "ampliconSize"),
    ("contaminationAssessment", MIQE + "contaminationAssessment"),
    ("amplificationReactionVolume", MIQE + "amplificationReactionVolume"),
    ("amplificationReactionVolumeUnit", MIQE + "amplificationReactionVolumeUnit"),
    ("pcr_analysis_software", MIQE + "pcr_analysis_software"),
])

# FAIRe checklist terms (FAIR-eDNA/FAIRe_checklist; Takahashi et al. 2025, Environmental DNA
# 7:e70100). Those that are also Darwin Core or DNA-derived extension terms go to their
# columns; the others are written to dynamicProperties under their FAIRe names.
FAIRE_TERMS = frozenset("""
accessRights accession_id accession_id_ref_db adapter_forward adapter_reverse alt amp_vis_method
ampliconSize amplificationReactionVolume annealingTemp assay_name assay_type assay_validation
associatedSequences associated_resource automaticBaselineValue
automaticThresholdQuantificationCycle barcoding_pcr_appr baselineValue bibliographicCitation
bioinfo_method_additional biological_rep biological_rep_relation block_ref block_seq block_taxa
checkls_ver checksum_filename checksum_filename2 checksum_method chimera_check_method
chimera_check_param chlorophyll code_repo commercial_mm concentration concentration_method
concentration_unit confidence_score custom_mm dataGeneralizations date_ext demux_max_mismatch
demux_tool detected_notDetected detection_criteria diss_inorg_carb diss_inorg_carb_unit
diss_inorg_nitro diss_inorg_nitro_unit diss_org_carb diss_org_carb_unit diss_org_nitro
diss_org_nitro_unit diss_oxygen diss_oxygen_unit dna_cleanup_0_1 dna_cleanup_method dna_sequence
dna_store_loc efficiency elev env_broad_scale env_local_scale env_medium error_rate_cutoff
error_rate_tool error_rate_type estimatedNumberOfCopies estimatedNumberOfCopies_error
estimatedNumberOfCopies_error_type estimatedNumberOfCopies_method estimatedNumberOfCopies_unit
eventDurationValue filename filename2 filter_diameter filter_material filter_name
filter_passive_active_0_1 filter_surface_area geo_loc_name habitat_natural_artificial_0_1
host_height host_height_unit host_length host_length_unit host_life_stage host_species
host_tot_mass host_tot_mass_unit humidity inhibition_check inhibition_check_0_1
input_read_count institution instrument intercept lib_conc lib_conc_meth lib_conc_unit lib_id
lib_layout lib_screen license light_intensity lod_method loq_method merge_min_overlap merge_tool
mid_forward mid_reverse min_len_cutoff min_len_tool min_reads_cutoff min_reads_cutoff_unit
min_reads_tool mod_date neg_cont_0_1 neg_cont_type nitrate nitrate_unit nitrite nitrite_unit
nitro nitro_unit nucl_acid_amp nucl_acid_ext nucl_acid_ext_kit nucl_acid_ext_lysis
nucl_acid_ext_method_additional nucl_acid_ext_modify nucl_acid_ext_sep num_detect org_carb
org_carb_unit org_matter org_matter_unit org_nitro org_nitro_unit otu_clust_cutoff otu_clust_tool
otu_db otu_db_custom otu_final_description otu_num_tax_assigned otu_raw_description
otu_seq_comp_appr output_otu_num output_read_count parent_project_id part_org_carb
part_org_carb_unit part_org_nitro part_org_nitro_unit pcr2_amplificationReactionVolume
pcr2_analysis_software pcr2_annealingTemp pcr2_commercial_mm pcr2_cond pcr2_custom_mm pcr2_cycles
pcr2_dna_vol pcr2_method_additional pcr2_plate_id pcr2_thermocycler pcr_0_1 pcr_analysis_software
pcr_assay_lod pcr_assay_lod_LL pcr_assay_lod_UL pcr_assay_lod_techreps pcr_assay_lod_unit
pcr_assay_loq pcr_assay_loq_LL pcr_assay_loq_UL pcr_assay_loq_techreps pcr_assay_loq_unit
pcr_cond pcr_cycles pcr_dna_vol pcr_method_additional pcr_plate_id pcr_primer_conc_forward
pcr_primer_conc_reverse pcr_primer_forward pcr_primer_name_forward pcr_primer_name_reverse
pcr_primer_reference_forward pcr_primer_reference_reverse pcr_primer_reverse
pcr_primer_vol_forward pcr_primer_vol_reverse pcr_rep percent_match percent_query_cover ph
ph_meth phix_perc platform pool_dna_num pos_cont_0_1 pos_cont_type precip_chem_prep
precip_force_prep precip_temp_prep precip_time_prep prefilter_material prep_method_additional
prepped_samp_store_dur prepped_samp_store_sol prepped_samp_store_temp probeQuencher probeReporter
probe_conc probe_ref probe_seq project_contact project_id project_name pump_flow_rate
pump_flow_rate_unit quantificationCycle r2 ratioOfAbsorbance260_280 rel_cont_id rfu salinity
samp_category samp_collect_device samp_collect_method samp_mat_process samp_name samp_size
samp_size_unit samp_store_dur samp_store_loc samp_store_method_additional samp_store_sol
samp_store_temp samp_vol_we_dna_ext samp_vol_we_dna_ext_unit samp_weather sample_composed_of
sample_derived_from screen_contam_method screen_geograph_method screen_nontarget_method
screen_other seq_id seq_kit seq_method_additional seq_run_id sequencing_location size_frac
size_frac_low slope solar_irradiance sop_bioinformatics stationed_sample_dur std_conc
std_conc_unit std_seq std_source std_type sterilise_method study_factor suspend_part_matter
targetTaxonomicAssay targetTaxonomicScope target_gene target_subfragment
targeted_detection_method_additional tax_assign_cat tax_class_collapse tax_class_id_cutoff
tax_class_other tax_class_query_cutoff taxonID_db technical_rep_id temp thermocycler
thresholdQuantificationCycle tidal_stage tot_carb tot_carb_unit tot_depth_water_col
tot_diss_nitro tot_diss_nitro_unit tot_inorg_nitro tot_inorg_nitro_unit tot_nitro
tot_nitro_cont_meth tot_nitro_content tot_nitro_content_unit tot_nitro_unit tot_org_c_meth
tot_org_carb tot_org_carb_unit tot_part_carb tot_part_carb_unit trim_method trim_param turbidity
water_current wind_direction wind_speed
""".split())
# FAIRe terms already carried by another column of the archive (not repeated in dynamicProperties)
FAIRE_MAPPED = {"samp_collect_method": "samp_collec_method", "samp_collect_device": "samp_collec_device",
                "geo_loc_name": "locality", "dna_sequence": "DNA_sequence"}
# FAIRe terms written to dynamicProperties: everything without a DwC / DNA-extension column
DYN_TERMS = sorted(t for t in FAIRE_TERMS
                   if t not in DNA_TERMS and t not in OCCURRENCE_TERMS and t not in FAIRE_MAPPED)
# columns of the OTU table that hold per-OTU values instead of read counts
OTU_TERM_COLUMNS = set(DNA_TERMS) | set(OCCURRENCE_TERMS) | FAIRE_TERMS
DNA_ROWTYPE = "http://rs.gbif.org/terms/1.0/DNADerivedData"

# Fields of the DNA extension that GBIF/OBIS treat as essential / strongly recommended
DNA_ESSENTIAL = ["DNA_sequence", "target_gene", "pcr_primer_forward", "pcr_primer_reverse",
                 "pcr_primer_name_forward", "pcr_primer_name_reverse", "pcr_primer_reference",
                 "env_medium", "seq_meth", "otu_db", "otu_seq_comp_appr", "sop"]

# NCBI rank -> Darwin Core column. These are all the rank terms Darwin Core defines (and that
# GBIF and OBIS index). Intermediate NCBI ranks without a DwC term (subphylum, superclass,
# subclass, infraclass, superorder, suborder, infraorder, cohort, clades, ...) are kept in
# higherClassification only.
RANK_MAP = {"kingdom": "kingdom", "phylum": "phylum", "class": "class", "order": "order",
            "superfamily": "superfamily", "family": "family", "subfamily": "subfamily",
            "tribe": "tribe", "subtribe": "subtribe", "genus": "genus", "subgenus": "subgenus",
            "species": "species", "subspecies": "subspecies", "varietas": "subspecies",
            "forma": "subspecies"}
TAXON_COLUMNS = ["scientificName", "taxonRank", "taxonID", "identificationQualifier",
                 "kingdom", "phylum", "class", "order", "superfamily", "family", "subfamily",
                 "tribe", "subtribe", "genus", "genericName", "subgenus", "infragenericEpithet",
                 "specificEpithet", "infraspecificEpithet", "higherClassification"]
# Open-nomenclature / placeholder species names (e.g. 'Pegusa sp. XYZ', 'uncultured Rissoa')
OPEN_NOMEN = re.compile(r"\b(sp|spp|cf|aff|nr|gr)\.|\b(environmental|uncultured|unidentified|"
                        r"unclassified|metagenome)\b", re.I)
# NCBI kingdom names -> names matched by the GBIF backbone / WoRMS
KINGDOM_FIX = {"Metazoa": "Animalia", "Viridiplantae": "Plantae"}
# NCBI "domain" used as fallback kingdom for prokaryotes
DOMAIN_AS_KINGDOM = {"Bacteria": "Bacteria", "Archaea": "Archaea"}

# Rough bounding boxes (lon_min, lat_min, lon_max, lat_max) used only for a plausibility check
WATER_BODY_BBOX = {
    "north sea": (-4.5, 50.9, 12.0, 62.0),
    "baltic sea": (9.0, 53.0, 30.5, 66.0),
    "bay of biscay": (-10.0, 43.2, -1.0, 48.5),
    "english channel": (-6.0, 48.4, 2.0, 51.2),
    "mediterranean sea": (-6.0, 30.0, 36.5, 46.0),
    "aegean sea": (22.0, 35.0, 28.5, 41.0),
    "adriatic sea": (12.0, 39.8, 20.0, 45.9),
    "black sea": (27.0, 40.8, 42.0, 47.0),
    "irish sea": (-6.5, 52.0, -2.8, 55.0),
    "celtic sea": (-11.0, 48.0, -5.0, 52.5),
}

FISH_CLADES = {"Actinopterygii", "Chondrichthyes", "Dipnoi", "Coelacanthiformes", "Cyclostomata",
               "Myxini", "Hyperoartia", "Actinopteri", "Elasmobranchii", "Holocephali"}

CONTROL_WORDS = ("control", "blank", "standard", "mock", "neg", "pos")


# --------------------------------------------------------------------------------------
# QC report collector
# --------------------------------------------------------------------------------------
class Report:
    def __init__(self):
        self.items = {"ERROR": [], "WARNING": [], "INFERRED": [], "INFO": []}

    def add(self, level, msg):
        if msg not in self.items[level]:
            self.items[level].append(msg)

    def markdown(self, header):
        out = [header, ""]
        titles = {"ERROR": "Errors (fix before publishing)",
                  "WARNING": "Warnings (missing or suspicious information)",
                  "INFERRED": "Values inferred automatically (please verify)",
                  "INFO": "Summary"}
        for lvl in ("ERROR", "WARNING", "INFERRED", "INFO"):
            out.append(f"## {titles[lvl]} ({len(self.items[lvl])})")
            out.append("")
            out += [f"- {m}" for m in self.items[lvl]] or ["- none"]
            out.append("")
        return "\n".join(out)


# --------------------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------------------
def sniff_read(path):
    """Read CSV/TSV with BOM and CRLF handling; returns list of dicts and fieldnames.
    Comma, tab and semicolon separators are recognised. Rows whose first cell starts with '#'
    (guidance rows of the metadata template builder) are skipped."""
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        text = fh.read()
    first = text.splitlines()[0] if text else ""
    delim = max(("\t", ",", ";"), key=first.count)
    reader = csv.DictReader(io.StringIO(text), delimiter=delim)
    fields = [f.strip() for f in reader.fieldnames]
    rows = []
    for r in reader:
        row = {k.strip(): (v.strip() if isinstance(v, str) else "") for k, v in r.items()
               if k is not None}
        if fields and row.get(fields[0], "").startswith("#"):
            continue
        rows.append(row)
    return rows, fields


def read_fasta(path):
    seqs, cur = {}, None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                cur = re.split(r"[\s;]", line[1:], maxsplit=1)[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


# --------------------------------------------------------------------------------------
# Parsers / helpers
# --------------------------------------------------------------------------------------
def parse_lineage(s):
    """Parse 'rank_Name_taxid;rank_Name_taxid;...' (NCBI-style, as written by e.g. OPTICS/vsearch)."""
    s = (s or "").strip()
    if not s or s.lower() in ("unassigned", "na", "none", "no_hit", "nohit", "unclassified"):
        return None
    nodes = []
    for tok in s.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"^([A-Za-z]+)_(.+?)(?:_(\d+))?$", tok)
        if not m:
            nodes.append({"rank": "unknown", "name": tok, "taxid": ""})
            continue
        rank, name, taxid = m.group(1).lower(), m.group(2), m.group(3) or ""
        if rank == "root" and name.isdigit():
            continue
        if rank in ("cellularroot",):
            continue
        nodes.append({"rank": rank, "name": name.replace("_", " ").strip(), "taxid": taxid})
    return nodes or None


def lineage_to_dwc(nodes):
    """Split a parsed lineage into the Darwin Core taxonomy columns used by GBIF and OBIS."""
    out = {c: "" for c in TAXON_COLUMNS}
    by_rank = {}
    for n in nodes:
        col = RANK_MAP.get(n["rank"])
        if col:
            by_rank[col] = n
        if n["rank"] in ("domain", "superkingdom") and "kingdom" not in by_rank \
                and n["name"] in DOMAIN_AS_KINGDOM:
            by_rank["kingdom"] = dict(n, name=DOMAIN_AS_KINGDOM[n["name"]])

    # Placeholder species names ('X sp. ABC', 'uncultured X') are not names GBIF/OBIS can match:
    # identify to the next ranked taxon up and record the qualifier.
    qualifier = ""
    for col in ("subspecies", "species"):
        n = by_rank.get(col)
        if n and OPEN_NOMEN.search(n["name"]):
            m = OPEN_NOMEN.search(n["name"])
            qualifier = (m.group(1).lower() + ".") if m.group(1) else m.group(2).lower()
            del by_rank[col]

    for col in ("kingdom", "phylum", "class", "order", "superfamily", "family",
                "subfamily", "tribe", "subtribe", "genus"):
        if col in by_rank:
            out[col] = by_rank[col]["name"]
    if out["kingdom"]:
        out["kingdom"] = KINGDOM_FIX.get(out["kingdom"], out["kingdom"])

    sp = by_rank.get("species", {}).get("name", "")
    ssp = by_rank.get("subspecies", {}).get("name", "")
    if sp and " " in sp:
        g, ep = sp.split(" ", 1)
        out["genus"] = out["genus"] or g
        out["specificEpithet"] = ep.split(" ")[0]
    if ssp:
        parts = re.sub(r"\b(subsp|ssp|var|f)\.\s*", "", ssp).split()
        if len(parts) >= 3:
            out["infraspecificEpithet"] = parts[-1]
            out["genus"] = out["genus"] or parts[0]
            out["specificEpithet"] = out["specificEpithet"] or parts[1]
    if "subgenus" in by_rank:
        sg = re.sub(r"^.*\((.+)\)$", r"\1", by_rank["subgenus"]["name"]).split()[-1]
        out["infragenericEpithet"] = sg
        out["subgenus"] = f"{out['genus']} ({sg})" if out["genus"] else sg

    # lowest DwC-ranked node becomes the identification
    order = ["kingdom", "phylum", "class", "order", "superfamily", "family", "subfamily",
             "tribe", "subtribe", "genus", "subgenus", "species", "subspecies"]
    lowest = next((c for c in reversed(order) if c in by_rank), None)
    if lowest:
        leaf = by_rank[lowest]
        out["taxonRank"] = leaf["rank"]
    else:  # nothing at a DwC rank: fall back to the deepest named node
        ranked = [n for n in nodes if n["rank"] not in ("clade", "norank", "unknown", "root")
                  and "incertae sedis" not in n["name"].lower()]
        leaf = ranked[-1] if ranked else nodes[-1]
        out["taxonRank"] = leaf["rank"] if leaf["rank"] not in ("clade", "norank", "unknown") else ""
    out["scientificName"] = out["subgenus"] if lowest == "subgenus" else leaf["name"]
    out["taxonID"] = f"NCBI:txid{leaf['taxid']}" if leaf["taxid"] else ""
    out["identificationQualifier"] = qualifier
    if lowest in ("genus", "subgenus", "species", "subspecies"):
        out["genericName"] = out["genus"]
    out["higherClassification"] = " | ".join(n["name"] for n in nodes[:nodes.index(leaf)])
    return out


DATE_FORMATS = [("%Y-%m-%d", False), ("%Y-%m-%dT%H:%M:%S", False), ("%Y-%m-%dT%H:%M", False),
                ("%d.%m.%Y", True), ("%d/%m/%Y", True), ("%d-%m-%Y", True),
                ("%Y/%m/%d", False), ("%Y%m%d", False)]


def parse_date(raw, report):
    raw = (raw or "").strip()
    if not raw:
        return ""
    if re.match(r"^\d{4}(-\d{2}(-\d{2})?)?(/\d{4}(-\d{2}(-\d{2})?)?)?$", raw):
        return raw  # already ISO 8601 (incl. ranges / partial dates)
    for fmt, day_first in DATE_FORMATS:
        try:
            d = dt.datetime.strptime(raw, fmt)
        except ValueError:
            continue
        iso = d.strftime("%Y-%m-%dT%H:%M:%S") if "H" in fmt else d.date().isoformat()
        msg = f"eventDate '{raw}' converted to ISO 8601 '{iso}'"
        if day_first and d.day <= 12 and d.month <= 12 and d.day != d.month:
            msg += " (assumed DD.MM.YYYY - ambiguous, please confirm)"
        report.add("INFERRED", msg)
        return iso
    report.add("ERROR", f"eventDate '{raw}' could not be parsed; supply ISO 8601 (YYYY-MM-DD)")
    return raw


def to_float(v):
    try:
        return float(str(v).replace(",", "."))
    except (TypeError, ValueError):
        return None


def is_primer_seq(s):
    return bool(s) and re.fullmatch(r"[ACGTURYSWKMBDHVNI]+", s.upper().replace(" ", "")) is not None


def infer_target_gene(assay_name, primer_names):
    blob = " ".join([assay_name or ""] + list(primer_names)).upper()
    for key, gene in [("12S", "12S rRNA"), ("16S", "16S rRNA"), ("18S", "18S rRNA"),
                      ("28S", "28S rRNA"), ("COI", "COI"), ("CO1", "COI"), ("ITS", "ITS"),
                      ("RBCL", "rbcL"), ("CYTB", "cytb"), ("MATK", "matK")]:
        if key in blob:
            return gene
    return ""


def is_fish(nodes):
    return any(n["name"] in FISH_CLADES for n in nodes)


# --------------------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------------------
def write_tsv(rows, fields):
    clean = lambda v: re.sub(r"[\t\r\n]+", " ", "" if v is None else str(v))  # noqa: E731
    lines = ["\t".join(fields)]
    lines += ["\t".join(clean(r.get(f, "")) for f in fields) for r in rows]
    return "\n".join(lines) + "\n"


def meta_xml(occ_fields, dna_fields):
    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<archive xmlns="http://rs.tdwg.org/dwc/text/" metadata="eml.xml">',
             f'  <core encoding="UTF-8" fieldsTerminatedBy="\\t" linesTerminatedBy="\\n" '
             f'fieldsEnclosedBy="" ignoreHeaderLines="1" rowType="{OCC_URI["rowType"]}">',
             '    <files><location>occurrence.txt</location></files>',
             '    <id index="0"/>']
    for i, f in enumerate(occ_fields):
        lines.append(f'    <field index="{i}" term="{OCCURRENCE_TERMS[f]}{f}"/>')
    lines.append('  </core>')
    if dna_fields:
        lines += [f'  <extension encoding="UTF-8" fieldsTerminatedBy="\\t" linesTerminatedBy="\\n" '
                  f'fieldsEnclosedBy="" ignoreHeaderLines="1" rowType="{DNA_ROWTYPE}">',
                  '    <files><location>dna_derived_data.txt</location></files>',
                  '    <coreid index="0"/>']
        for i, f in enumerate(dna_fields[1:], start=1):
            lines.append(f'    <field index="{i}" term="{DNA_TERMS[f]}"/>')
        lines.append('  </extension>')
    lines.append('</archive>')
    return "\n".join(lines) + "\n"


def eml_xml(cfg, cov, methods_text):
    ds = cfg.get("dataset", {})
    e = lambda s: escape(str(s or ""))  # noqa: E731

    def party(tag, p):
        ind = ""
        if p.get("givenName") or p.get("surName"):
            ind = (f"<individualName><givenName>{e(p.get('givenName'))}</givenName>"
                   f"<surName>{e(p.get('surName') or p.get('givenName'))}</surName></individualName>")
        org = f"<organizationName>{e(p.get('organization'))}</organizationName>" if p.get("organization") else ""
        mail = f"<electronicMailAddress>{e(p.get('email'))}</electronicMailAddress>" if p.get("email") else ""
        orcid = (f'<userId directory="https://orcid.org/">{e(p.get("orcid"))}</userId>'
                 if p.get("orcid") else "")
        return f"<{tag}>{ind}{org}{mail}{orcid}</{tag}>"

    creators = cfg.get("creators") or [{"givenName": "", "surName": "UNKNOWN"}]
    contact = cfg.get("contact") or creators[0]
    kw = "".join(f"<keyword>{e(k)}</keyword>" for k in ds.get("keywords", []))
    lic = ds.get("license", "CC-BY 4.0")
    lic_url = {"CC0 1.0": "http://creativecommons.org/publicdomain/zero/1.0/legalcode",
               "CC-BY 4.0": "http://creativecommons.org/licenses/by/4.0/legalcode",
               "CC-BY-NC 4.0": "http://creativecommons.org/licenses/by-nc/4.0/legalcode"}.get(lic, "")
    rights = (f'This work is licensed under a <ulink url="{lic_url}"><citetitle>{e(lic)}</citetitle>'
              f'</ulink>.' if lic_url else e(lic))
    tax = "".join(f"<taxonomicClassification><taxonRankName>{e(r)}</taxonRankName>"
                  f"<taxonRankValue>{e(v)}</taxonRankValue></taxonomicClassification>"
                  for r, v in cov["taxa"])
    temporal = ""
    if cov["dates"]:
        d0, d1 = min(cov["dates"])[:10], max(cov["dates"])[:10]
        temporal = (f"<temporalCoverage><rangeOfDates><beginDate><calendarDate>{d0}</calendarDate>"
                    f"</beginDate><endDate><calendarDate>{d1}</calendarDate></endDate>"
                    f"</rangeOfDates></temporalCoverage>")
    geo = ""
    if cov["lats"]:
        geo = (f"<geographicCoverage><geographicDescription>{e(', '.join(sorted(cov['places'])) or 'see records')}"
               f"</geographicDescription><boundingCoordinates>"
               f"<westBoundingCoordinate>{min(cov['lons'])}</westBoundingCoordinate>"
               f"<eastBoundingCoordinate>{max(cov['lons'])}</eastBoundingCoordinate>"
               f"<northBoundingCoordinate>{max(cov['lats'])}</northBoundingCoordinate>"
               f"<southBoundingCoordinate>{min(cov['lats'])}</southBoundingCoordinate>"
               f"</boundingCoordinates></geographicCoverage>")
    steps = "".join(f"<methodStep><description><para>{e(s)}</para></description></methodStep>"
                    for s in methods_text)
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<eml:eml xmlns:eml="eml://ecoinformatics.org/eml-2.1.1"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="eml://ecoinformatics.org/eml-2.1.1 http://rs.gbif.org/schema/eml-gbif-profile/1.2/eml.xsd"
         packageId="{e(ds.get('packageId') or uuid.uuid4())}" system="http://gbif.org" scope="system" xml:lang="{e(ds.get('language', 'eng'))}">
<dataset>
  <title xml:lang="{e(ds.get('language', 'eng'))}">{e(ds.get('title'))}</title>
  {''.join(party('creator', c) for c in creators)}
  {party('metadataProvider', cfg.get('metadataProvider') or creators[0])}
  <pubDate>{dt.date.today().isoformat()}</pubDate>
  <language>{e(ds.get('language', 'eng'))}</language>
  <abstract><para>{e(ds.get('abstract'))}</para></abstract>
  <keywordSet>{kw}<keywordThesaurus>N/A</keywordThesaurus></keywordSet>
  <intellectualRights><para>{rights}</para></intellectualRights>
  <coverage>{geo}{temporal}<taxonomicCoverage>{tax}</taxonomicCoverage></coverage>
  {party('contact', contact)}
  <methods>{steps}<sampling><studyExtent><description><para>{e(ds.get('studyExtent') or 'Environmental DNA metabarcoding survey')}</para></description></studyExtent><samplingDescription><para>{e(ds.get('samplingDescription') or 'See DNA-derived data extension (MIxS/FAIRe fields) for sample processing details.')}</para></samplingDescription></sampling></methods>
  <project><title>{e(ds.get('projectTitle') or cov.get('project') or ds.get('title'))}</title>{party('personnel', contact).replace('</personnel>', '<role>pointOfContact</role></personnel>')}</project>
</dataset>
</eml:eml>
"""


# --------------------------------------------------------------------------------------
# Main conversion
# --------------------------------------------------------------------------------------
CONFIG_TEMPLATE = {
    "dataset": {
        "title": "eDNA metabarcoding of <target> in <area>, <year>",
        "abstract": "Describe aim, sampling, marker, bioinformatics and reference database.",
        "license": "CC-BY 4.0",
        "language": "eng",
        "keywords": ["eDNA", "metabarcoding", "DNA-derived data"],
        "projectTitle": "",
        "datasetName": "",
        "institutionCode": "",
    },
    "creators": [{"givenName": "", "surName": "", "organization": "", "email": "", "orcid": ""}],
    "contact": None,
    "defaults": {
        "_comment": "Used when the sample metadata lacks the column. Metadata values always win.",
        "geodeticDatum": "EPSG:4326",
        "coordinateUncertaintyInMeters": "",
        "countryCode": "",
        "samplingProtocol": "",
        "target_gene": "",
        "target_subfragment": "",
        "pcr_primer_forward": "",
        "pcr_primer_reverse": "",
        "pcr_primer_reference": "",
        "seq_meth": "",
        "lib_layout": "",
        "nucl_acid_ext": "",
        "samp_vol_we_dna_ext": "",
        "size_frac": "",
        "otu_class_appr": "",
        "otu_seq_comp_appr": "",
        "otu_db": "",
        "sop": "",
        "_comment2": "Any other DNA-derived extension term (e.g. seq_quality_check, chimera_check, "
                     "ref_db, sim_search_meth, tax_class) or FAIRe term (e.g. trim_method, "
                     "otu_clust_tool, tax_class_id_cutoff) can be given here too.",
    },
    "sample_values": {"_comment": "per-sample values from the pipeline, used where the metadata has "
                                  "none: {\"sample1\": {\"lib_reads_seqd\": 12345, \"output_read_count\": 12000}}"},
    "primer_sequences": {"_comment": "primer name -> 5'-3' sequence", "16S_F": "", "16S_R": ""},
}


def deep_get(d, *keys, default=""):
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
    return default if d is None else d


def convert(args):
    rep = Report()
    cfg = {}
    if args.config:
        with open(args.config, encoding="utf-8") as fh:
            cfg = json.load(fh)
    defaults = {k: v for k, v in deep_get(cfg, "defaults", default={}).items() if not k.startswith("_")}
    primer_seqs = {k: v for k, v in deep_get(cfg, "primer_sequences", default={}).items()
                   if not k.startswith("_") and v}

    # ---------- read inputs ----------
    otu_rows, otu_fields = sniff_read(args.otu)
    meta_rows, meta_fields = sniff_read(args.metadata)
    seqs = read_fasta(args.fasta) if args.fasta else {}

    id_col = otu_fields[0]
    term_cols = [f for f in otu_fields[1:] if f in OTU_TERM_COLUMNS]   # per-OTU term values
    other = [f for f in otu_fields if f not in term_cols]
    tax_col = next((f for f in other if f.lower() in ("taxon", "taxonomy", "lineage",
                                                        "classification", "tax")), None)
    ident_col = next((f for f in other if "ident" in f.lower()), None)
    ref_col = next((f for f in other if "ref" in f.lower()), None)
    sample_cols = [f for f in other if f not in (id_col, tax_col, ident_col, ref_col)]
    if not tax_col:
        rep.add("ERROR", "No taxonomy column found in OTU table (expected 'taxon'/'taxonomy').")

    name_col = next((f for f in meta_fields if f.lower() in ("samp_name", "sample", "sample_id",
                                                               "materialsampleid")), meta_fields[0])
    meta = {r[name_col]: r for r in meta_rows}
    sample_values = {k: v for k, v in deep_get(cfg, "sample_values", default={}).items()
                     if not k.startswith("_") and isinstance(v, dict)}
    for s, r in meta.items():
        for k, v in sample_values.get(s, {}).items():     # pipeline values fill empty metadata
            if not str(r.get(k) or "").strip() and str(v).strip():
                r[k] = str(v)
        for faire, dwc in FAIRE_MAPPED.items():            # FAIRe name -> DNA-derived / MIxS name
            if dwc in DNA_TERMS and not r.get(dwc) and r.get(faire):
                r[dwc] = r[faire]

    # ---------- sample matching ----------
    in_table = set(sample_cols)
    in_meta = set(meta)
    for s in sorted(in_table - in_meta):
        rep.add("ERROR", f"Sample '{s}' is in the OTU table but has no metadata row - skipped.")
    for s in sorted(in_meta - in_table):
        rep.add("WARNING", f"Sample '{s}' has metadata but no column in the OTU table.")
    samples = [s for s in sample_cols if s in meta]

    # ---------- per-sample checks ----------
    totals = {s: 0 for s in samples}
    for r in otu_rows:
        for s in samples:
            v = to_float(r.get(s)) or 0
            totals[s] += int(v)

    event_ids = Counter(meta[s].get("eventID", "") for s in samples)
    controls = []
    sample_info = {}
    grouped = OrderedDict()   # message -> [samples]
    for s in samples:
        m = meta[s]
        info = {}
        cat = m.get("samp_category", "sample").strip().lower()
        if any(w in cat for w in CONTROL_WORDS) and cat != "sample":
            controls.append(s)
        lat, lon = to_float(m.get("decimalLatitude")), to_float(m.get("decimalLongitude"))
        if lat is None or lon is None:
            rep.add("ERROR", f"{s}: decimalLatitude/decimalLongitude missing or not numeric.")
        elif not (-90 <= lat <= 90 and -180 <= lon <= 180):
            rep.add("ERROR", f"{s}: coordinates out of range ({lat}, {lon}).")
        else:
            place = (m.get("geo_loc_name") or "").lower()
            for wb, (x0, y0, x1, y1) in WATER_BODY_BBOX.items():
                if wb in place and not (x0 <= lon <= x1 and y0 <= lat <= y1):
                    hit = [n for n, (a, b, c, d) in WATER_BODY_BBOX.items()
                           if a <= lon <= c and b <= lat <= d and n != "mediterranean sea"]
                    grouped.setdefault(f"coordinates ({lat}, {lon}) fall outside '{m.get('geo_loc_name')}'"
                            + (f" - they look like {', '.join(h.title() for h in hit)}" if hit else "")
                            + ". Check for swapped/mis-signed lat/lon or wrong place name.", []).append(s)
            if "belgium" in (m.get("eventID", "") + place).lower() and not (2.0 <= lon <= 6.5 and 49.4 <= lat <= 52.0):
                grouped.setdefault(f"eventID/place mention Belgium but coordinates ({lat}, {lon}) "
                                   "are not in/near Belgium.", []).append(s)
        info["eventDate"] = parse_date(m.get("eventDate", ""), rep)
        if not info["eventDate"]:
            rep.add("ERROR", f"{s}: eventDate missing.")
        if totals[s] == 0:
            rep.add("WARNING", f"{s}: zero reads in OTU table.")
        sample_info[s] = info

    for msg, ss in grouped.items():
        rep.add("WARNING", f"{', '.join(ss)}: {msg}")
    for eid, n in event_ids.items():
        if eid and n > 1:
            rep.add("INFO", f"eventID '{eid}' is shared by {n} samples - treated as replicate samples "
                    "of one sampling event (materialSampleID distinguishes them). If they are different "
                    "stations/times, give each its own eventID.")
        if not eid:
            rep.add("WARNING", "Some samples lack an eventID; one is generated from project_id + samp_name.")
    if controls:
        rep.add("INFO", f"Control samples excluded from occurrences: {', '.join(controls)} "
                "(their reads should be used for decontamination before publishing).")
    else:
        rep.add("WARNING", "No negative/positive controls listed in metadata (samp_category). "
                "eDNAqua-Plan and FAIRe recommend reporting field, extraction and PCR blanks.")

    # ---------- assay-level fields (FAIRe -> MIxS) ----------
    def mval(s, key):
        return str(meta[s].get(key) or defaults.get(key) or "").strip()

    first = meta[samples[0]] if samples else {}
    for p_key in ("pcr_primer_forward", "pcr_primer_reverse"):
        v = first.get(p_key, "") or defaults.get(p_key, "")
        if v and not is_primer_seq(v):
            if v in primer_seqs:
                rep.add("INFERRED", f"{p_key}: name '{v}' resolved to sequence via config primer_sequences.")
            else:
                rep.add("WARNING", f"{p_key} = '{v}' is a primer name, not a 5'-3' sequence. It is written to "
                        f"pcr_primer_name_{p_key.split('_')[-1]}; add the sequence via config 'primer_sequences'.")
    tg_existing = first.get("target_gene") or defaults.get("target_gene")
    tg_inferred = "" if tg_existing else infer_target_gene(
        first.get("assay_name", ""), [first.get("pcr_primer_forward", ""), first.get("pcr_primer_reverse", "")])
    if tg_inferred:
        rep.add("INFERRED", f"target_gene '{tg_inferred}' inferred from assay/primer names.")
    elif not tg_existing:
        rep.add("WARNING", "target_gene unknown - add it to metadata or config defaults.")

    if not seqs:
        rep.add("WARNING", "No OTU/ASV sequences supplied (--fasta). DNA_sequence is the key field of the "
                "DNA-derived extension; without it records cannot be re-identified later and "
                "unassigned OTUs carry almost no information. Export the representative sequences "
                "from the pipeline and rerun.")
    for k in DNA_ESSENTIAL:
        if k in ("DNA_sequence", "pcr_primer_forward", "pcr_primer_reverse",
                 "pcr_primer_name_forward", "pcr_primer_name_reverse", "target_gene"):
            continue
        if not (first.get(k) or defaults.get(k)):
            rep.add("WARNING", f"Recommended DNA-derived field '{k}' is empty (add to metadata or config defaults).")
    for k in ("coordinateUncertaintyInMeters", "countryCode", "samplingProtocol",
              "samp_vol_we_dna_ext", "size_frac", "nucl_acid_ext"):
        if not (first.get(k) or defaults.get(k)):
            rep.add("WARNING", f"Recommended field '{k}' is empty.")
    if not (first.get("geodeticDatum") or args.config and defaults.get("geodeticDatum")):
        rep.add("INFERRED", "geodeticDatum assumed EPSG:4326 (WGS84).")
    for k in ("title", "abstract"):
        if not deep_get(cfg, "dataset", k) or "<" in deep_get(cfg, "dataset", k):
            rep.add("WARNING", f"EML dataset {k} missing or still a template - supply --config.")
    if not deep_get(cfg, "creators"):
        rep.add("WARNING", "No creators in config; EML creator is a placeholder (GBIF/OBIS require one).")

    target = (first.get("targetTaxonomicAssay") or "").lower()

    # ---------- build occurrences ----------
    occ, dna = [], []
    tax_cache, off_target, unassigned_reads, n_unassigned = {}, set(), 0, 0
    project = first.get("project_id", "") or first.get("project_name", "")
    for r in otu_rows:
        otu = r[id_col]
        lineage = parse_lineage(r.get(tax_col, "")) if tax_col else None
        if lineage:
            tx = tax_cache.setdefault(otu, lineage_to_dwc(lineage))
            if "fish" in target and not is_fish(lineage):
                off_target.add(f"{otu} ({tx['scientificName']})")
        else:
            tx = {"scientificName": args.unassigned_name, "taxonRank": "", "taxonID": "",
                  "kingdom": "incertae sedis" if args.unassigned_name.lower() == "incertae sedis" else ""}
        ident = to_float(r.get(ident_col)) if ident_col else None
        refs = (r.get(ref_col) or "").strip() if ref_col else ""
        if lineage:
            remarks = f"Assigned with {args.classifier}"
            if ident:
                remarks += f"; {getattr(args, 'identity_label', '') or 'mean identity to references'} {ident:.3f}"
            if refs:
                remarks += f"; best reference hits: {refs.replace(';', ', ')}"
        else:
            remarks = f"No taxonomic assignment ({args.classifier}); OTU {otu}"
        otu_vals = {k: str(r.get(k) or "").strip() for k in term_cols}
        otu_vals = {k: v for k, v in otu_vals.items() if v and v.lower() != "nan"}

        for s in samples:
            reads = int(to_float(r.get(s)) or 0)
            if reads < args.min_reads or s in controls:
                continue
            if not lineage:
                unassigned_reads += reads
                if args.unassigned == "drop":
                    continue
                n_unassigned += 1
            m = meta[s]
            event_id = m.get("eventID") or f"{project}_{s}"
            assay = m.get("assay_name", "") or "assay"
            prefix = assay if project and assay.startswith(project) else f"{project}_{assay}"
            occ_id = f"{prefix}_{s}_{otu}".replace(" ", "_")
            place = m.get("geo_loc_name", "")
            dyn = {k: mval(s, k) for k in DYN_TERMS if mval(s, k)}       # FAIRe terms
            dyn.update({k: v for k, v in otu_vals.items() if k in DYN_TERMS})
            dyn["seq_id"] = otu
            row = {
                "occurrenceID": occ_id,
                "basisOfRecord": "MaterialSample",
                "occurrenceStatus": "present",
                "eventID": event_id,
                "parentEventID": m.get("parentEventID", ""),
                "eventDate": sample_info[s]["eventDate"],
                "verbatimEventDate": m.get("eventDate", ""),
                "materialSampleID": s,
                "recordedBy": m.get("recordedBy", ""),
                "samplingProtocol": mval(s, "samplingProtocol"),
                "locality": place,
                "waterBody": place if re.search(r"\b(sea|ocean|gulf|bay|strait|channel|lake|river|lagoon)\b",
                                                place, re.I) else "",
                "countryCode": mval(s, "countryCode"),
                "decimalLatitude": m.get("decimalLatitude", ""),
                "decimalLongitude": m.get("decimalLongitude", ""),
                "geodeticDatum": mval(s, "geodeticDatum") or "EPSG:4326",
                "coordinateUncertaintyInMeters": mval(s, "coordinateUncertaintyInMeters"),
                "minimumDepthInMeters": m.get("minimumDepthInMeters", "") or m.get("depth", ""),
                "maximumDepthInMeters": m.get("maximumDepthInMeters", "") or m.get("depth", ""),
                "nameAccordingTo": "NCBI Taxonomy" if lineage else "",
                "verbatimIdentification": (r.get(tax_col) or "").strip() if lineage else "unassigned",
                "identificationRemarks": remarks,
                "identificationReferences": args.identification_reference or "",
                "organismQuantity": reads,
                "organismQuantityType": "DNA sequence reads",
                "sampleSizeValue": totals[s],
                "sampleSizeUnit": "DNA sequence reads",
                "datasetName": deep_get(cfg, "dataset", "datasetName") or deep_get(cfg, "dataset", "title"),
                "institutionCode": deep_get(cfg, "dataset", "institutionCode"),
                "associatedSequences": mval(s, "associatedSequences"),
                "dynamicProperties": json.dumps(dyn, ensure_ascii=False),
            }
            row.update(tx)
            for k, v in otu_vals.items():          # per-OTU core values (e.g. scientificNameID)
                if k in OCCURRENCE_TERMS and k not in ("occurrenceID", "dynamicProperties"):
                    row[k] = v
            occ.append(row)

            fwd = m.get("pcr_primer_forward", "") or defaults.get("pcr_primer_forward", "")
            rev = m.get("pcr_primer_reverse", "") or defaults.get("pcr_primer_reverse", "")
            drow = {
                "occurrenceID": occ_id,
                "samp_name": s,
                "project_name": project,
                "env_broad_scale": m.get("env_broad_scale", ""),
                "env_local_scale": m.get("env_local_scale", ""),
                "env_medium": m.get("env_medium", ""),
                "target_gene": mval(s, "target_gene") or tg_inferred,
                "pcr_primer_forward": fwd if is_primer_seq(fwd) else primer_seqs.get(fwd, ""),
                "pcr_primer_reverse": rev if is_primer_seq(rev) else primer_seqs.get(rev, ""),
                "pcr_primer_name_forward": m.get("pcr_primer_name_forward", "") or ("" if is_primer_seq(fwd) else fwd),
                "pcr_primer_name_reverse": m.get("pcr_primer_name_reverse", "") or ("" if is_primer_seq(rev) else rev),
                "DNA_sequence": seqs.get(otu, ""),
            }
            for k in DNA_TERMS:                     # every other extension term with a value
                if k not in drow:
                    drow[k] = otu_vals.get(k) or mval(s, k)
            dna.append(drow)

    if seqs:
        missing = {r[id_col] for r in otu_rows} - set(seqs)
        if missing:
            rep.add("WARNING", f"{len(missing)} OTUs have no sequence in the FASTA (e.g. {sorted(missing)[:5]}).")
    if off_target:
        rep.add("WARNING", f"Target taxon is '{target}' but {len(off_target)} OTU(s) are not fish: "
                f"{'; '.join(sorted(off_target))}. Kept as genuine detections - check for contamination "
                "or non-specific amplification.")
    low = [f"{o} ({tax_cache[o]['scientificName']}, {to_float(r.get(ident_col)):.3f})"
           for r in otu_rows for o in [r[id_col]]
           if o in tax_cache and ident_col and (to_float(r.get(ident_col)) or 1) < args.species_identity]
    if low:
        rep.add("WARNING", f"Assignments below {args.species_identity:.0%} identity (rank may be over-precise for "
                f"species/genus-level claims): {'; '.join(low)}. Consider raising them to a higher rank.")
    gaps = []
    for o, tx in tax_cache.items():
        want = ["kingdom", "phylum", "class", "order", "family", "genus"]
        depth = {"kingdom": 1, "phylum": 2, "class": 3, "order": 4, "superfamily": 4, "family": 5,
                 "subfamily": 5, "tribe": 5, "subtribe": 5, "genus": 6, "subgenus": 6,
                 "species": 6, "subspecies": 6}.get(tx["taxonRank"], 0)
        miss = [c for c in want[:depth] if not tx.get(c)]
        if miss:
            gaps.append(f"{o} ({tx['scientificName']}): no {', '.join(miss)} in lineage")
        if tx.get("identificationQualifier"):
            rep.add("INFERRED", f"{o}: placeholder species name in lineage; identified as "
                    f"{tx['scientificName']} ({tx['taxonRank']}) with identificationQualifier "
                    f"'{tx['identificationQualifier']}' (full name kept in verbatimIdentification).")
    if gaps:
        rep.add("WARNING", "Major ranks absent from the NCBI lineage (left empty; GBIF/OBIS fill them from "
                "their own backbone when the name matches): " + "; ".join(gaps) + ".")
    if any(t.get("kingdom") == "Animalia" for t in tax_cache.values()):
        rep.add("INFERRED", "NCBI kingdom 'Metazoa' mapped to 'Animalia' for GBIF/OBIS name matching.")

    total_reads = sum(totals[s] for s in samples if s not in controls)
    n_otu = len(otu_rows)
    n_assigned = len(tax_cache)
    rep.add("INFO", f"OTUs: {n_otu} ({n_assigned} assigned, {n_otu - n_assigned} unassigned).")
    rep.add("INFO", f"Samples: {len(samples)} ({len(controls)} controls). Total reads (non-control): {total_reads:,}.")
    if total_reads:
        rep.add("INFO", f"Reads in unassigned OTUs: {unassigned_reads:,} ({unassigned_reads / total_reads:.1%}).")
    rep.add("INFO", f"Occurrence records written: {len(occ)} "
            f"({n_unassigned if args.unassigned == 'keep' else 0} unassigned; unassigned mode = '{args.unassigned}').")

    # ---------- columns: drop fully empty optional ones ----------
    keep_always = {"occurrenceID", "basisOfRecord", "eventDate", "scientificName", "decimalLatitude",
                   "decimalLongitude", "occurrenceStatus", "organismQuantity", "organismQuantityType",
                   "sampleSizeValue", "sampleSizeUnit", "verbatimIdentification"} | set(TAXON_COLUMNS)
    occ_fields = [f for f in OCCURRENCE_TERMS if f in keep_always or any(r.get(f) not in ("", None) for r in occ)]
    dna_cols = [f for f in DNA_TERMS if any(r.get(f) for r in dna)]
    dna_fields = ["occurrenceID"] + dna_cols

    # ---------- EML ----------
    cov = {"lats": [], "lons": [], "dates": [], "places": set(), "taxa": [], "project": project}
    for s in samples:
        la, lo = to_float(meta[s].get("decimalLatitude")), to_float(meta[s].get("decimalLongitude"))
        if la is not None and lo is not None:
            cov["lats"].append(la); cov["lons"].append(lo)
        if sample_info[s]["eventDate"]:
            cov["dates"].append(sample_info[s]["eventDate"])
        if meta[s].get("geo_loc_name"):
            cov["places"].add(meta[s]["geo_loc_name"])
    seen = set()
    for tx in tax_cache.values():
        for rk in ("phylum", "class", "family"):
            if tx.get(rk) and (rk, tx[rk]) not in seen:
                seen.add((rk, tx[rk])); cov["taxa"].append((rk, tx[rk]))
    methods = [
        f"Sample metadata: {args.metadata.split('/')[-1]}; OTU table: {args.otu.split('/')[-1]}.",
        f"Assay {first.get('assay_name', '')} ({tg_existing or tg_inferred or 'marker unknown'}), primers "
        f"{first.get('pcr_primer_forward', '')}/{first.get('pcr_primer_reverse', '')}; target taxa: {target or 'n/a'}.",
        f"Taxonomic assignment: {args.classifier} against NCBI-taxonomy annotated references"
        + (f" ({defaults.get('otu_db')})" if defaults.get("otu_db") else "") + ".",
        "Occurrences = OTU x sample with reads >= "
        f"{args.min_reads}; organismQuantity = reads, sampleSizeValue = total reads per sample. "
        "Converted with edna2dwca.py following GBIF/OBIS guidance for DNA-derived data and eDNAqua-Plan recommendations.",
    ]
    faire_ds = [(k, str(defaults[k])) for k in DYN_TERMS if str(defaults.get(k) or "").strip()]
    if faire_ds:
        methods.append("Bioinformatics (FAIRe terms; also in dynamicProperties of every record): "
                       + "; ".join(f"{k} = {v}" for k, v in faire_ds) + ".")
    dyn_keys = sorted({k for r in occ for k in json.loads(r["dynamicProperties"])})
    rep.add("INFO", f"FAIRe terms in dynamicProperties ({len(dyn_keys)}): {', '.join(dyn_keys)}.")
    eml = eml_xml(cfg, cov, methods)

    # ---------- write ----------
    out_zip = args.output if args.output.endswith(".zip") else args.output + ".zip"
    with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("occurrence.txt", write_tsv(occ, occ_fields))
        z.writestr("dna_derived_data.txt", write_tsv(dna, dna_fields))
        z.writestr("meta.xml", meta_xml(occ_fields, dna_fields))
        z.writestr("eml.xml", eml)
    report_md = rep.markdown(f"# QC report - {out_zip.split('/')[-1]}\n\nGenerated "
                             f"{dt.datetime.now().isoformat(timespec='seconds')}")
    rep_path = out_zip[:-4] + "_QC_report.md"
    with open(rep_path, "w", encoding="utf-8") as fh:
        fh.write(report_md)

    print(f"Archive : {out_zip}  ({len(occ)} occurrences, {len(occ_fields)} core fields, "
          f"{len(dna_cols)} DNA-derived fields)")
    print(f"Report  : {rep_path}  ({len(rep.items['ERROR'])} errors, {len(rep.items['WARNING'])} warnings)")
    return 1 if rep.items["ERROR"] and args.strict else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--otu", help="OTU/ASV table (TSV/CSV): ID, one column per sample, taxon, [identity], [ref ids]")
    ap.add_argument("--metadata", help="Sample metadata (FAIRe-style CSV/TSV, one row per sample)")
    ap.add_argument("--config", help="JSON with EML info, defaults and primer sequences")
    ap.add_argument("--fasta", help="Representative OTU/ASV sequences (headers = OTU IDs)")
    ap.add_argument("-o", "--output", default="dwca_edna", help="Output archive name (.zip added)")
    ap.add_argument("--unassigned", choices=["keep", "drop"], default="keep",
                    help="Keep unassigned OTUs as records (default) or drop them")
    ap.add_argument("--unassigned-name", default="Biota",
                    help="scientificName for unassigned OTUs ('Biota' matches WoRMS/OBIS; or 'incertae sedis')")
    ap.add_argument("--min-reads", type=int, default=1, help="Minimum reads per OTU x sample (default 1)")
    ap.add_argument("--classifier", default="VSEARCH global alignment (OPTICS pipeline)",
                    help="Text describing the taxonomic assignment method")
    ap.add_argument("--identification-reference", default="", help="DOI/URL of pipeline or reference DB")
    ap.add_argument("--identity-label", default="mean identity to references",
                    help="How the identity column is described in identificationRemarks")
    ap.add_argument("--species-identity", type=float, default=0.97,
                    help="Flag assignments below this identity (default 0.97)")
    ap.add_argument("--strict", action="store_true", help="Exit code 1 if the QC report has errors")
    ap.add_argument("--write-config-template", metavar="PATH", help="Write a starter config JSON and exit")
    args = ap.parse_args()

    if args.write_config_template:
        with open(args.write_config_template, "w", encoding="utf-8") as fh:
            json.dump(CONFIG_TEMPLATE, fh, indent=2)
        print(f"Config template written to {args.write_config_template}")
        return 0
    if not (args.otu and args.metadata):
        ap.error("--otu and --metadata are required")
    return convert(args)


if __name__ == "__main__":
    sys.exit(main())
