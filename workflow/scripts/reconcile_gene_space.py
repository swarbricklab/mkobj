#!/usr/bin/env python3
"""
workflow/scripts/reconcile_gene_space.py
Lift an AnnData object's gene space from its source GENCODE build (the Cell Ranger
reference, e.g. GENCODE v32 / refdata-gex-GRCh38-2020-A) to a target GENCODE build
(e.g. v48, the gene set CELLxGENE schema 7.1.0 requires).

Ported from the brca_major_atlas CELLxGENE curation steps 02_retired_genes.py and
03_gene_map.py (the crosswalk) plus the apply half of 08_var.py / 09_write_object.py.

Key design fact
---------------
create_anndata.py sets ``adata.var_names`` to ``features['gene_id']`` — so the var
index is ALREADY the unversioned source Ensembl ID (ENSG…, the 2020-A set), in
count-matrix order, and ``var['gene_name']`` holds the symbol. The crosswalk is
therefore keyed directly on ``adata.var_names``; the publication's 01_gene_tables.py
reconstruction / make.unique machinery is not needed.

Structure
---------
``compute_gene_crosswalk(...)`` — steps 1-6 — is a pure function of the source IDs
and the three references, independent of any AnnData. It returns one row per source
gene, in input order, with the decision (keep / rename / drop), the target ID and the
provenance columns. The snakemake entry point (step 7) then applies that crosswalk to
the object. The validation gate in the test harness calls ``compute_gene_crosswalk``
directly to prove the port reproduces the publication's result.

A renamed column keeps its counts unchanged — two columns are never summed, and the
counts are not re-quantified against the target build.
"""

import gzip
import hashlib
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Constants (config-overridable; defaults reproduce the publication) ────────────
# Each gene span must cover at least this fraction of the other.
DEFAULT_MIN_RECIPROCAL_OVERLAP = 0.5
# At least this fraction of the source gene's exon bases must lie in the candidate.
DEFAULT_MIN_EXON_OVERLAP = 0.5
# A clone-based name, such as AL365255.1, is the accession of the clone sequence that a
# gene lies on, followed by a dot and a number. GENCODE gave such names to genes without
# an HGNC symbol and gives no target-build gene one.
CLONE_BASED_NAME = r"[A-Z]{1,2}\d{5,6}\.\d+"
# The columns of stable_id_event, as homo_sapiens_core_<rel>.sql defines them. No header row.
EVENT_COLUMNS = [
    "old_stable_id",
    "old_version",
    "new_stable_id",
    "new_version",
    "mapping_session_id",
    "type",
    "score",
]
# Each condition in the order tested, with the reason recorded when it is the first to fail.
REASONS = {
    "has_successor": "no successor in target build",
    "not_already_a_column": "merged into a gene already in the matrix",
    "single_candidate": "several history successors",
    "not_shared": "successor shared with another retired identifier",
    "same_name": "name differs in target build",
    "same_biotype": "biotype differs in target build",
    "same_locus": "candidate fails locus check",
    "exons_in_candidate": "candidate fails exon check",
}


# ── Checksum verification (no network; mirrors 00_download_references.sh) ──────────
def _md5(path):
    """Return the hex md5 of a file's contents (matches ``md5sum < file``)."""
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _bsd_sum(path):
    """Return the BSD checksum of a file (matches ``sum < file``), via coreutils ``sum``."""
    with open(path, "rb") as handle:
        result = subprocess.run(["sum"], stdin=handle, capture_output=True, text=True, check=True)
    # ``sum`` prints "<checksum> <blocks>"; normalise the whitespace.
    return " ".join(result.stdout.split())


def verify_file(path, md5=None, bsd_sum=None):
    """Stop with a clear error when a reference's checksum differs from the pinned one."""
    if md5 is not None:
        observed = _md5(path)
        if observed != md5:
            sys.exit(f"FAIL md5 of {path}: {observed} against {md5} (pinned checksum)")
        print(f"ok   md5 of {path}: {observed}")
    if bsd_sum is not None:
        observed = _bsd_sum(path)
        if observed != str(bsd_sum):
            sys.exit(f"FAIL BSD checksum of {path}: {observed} against {bsd_sum} (pinned checksum)")
        print(f"ok   BSD checksum of {path}: {observed}")


# ── GTF readers (ported verbatim from 01_gene_tables.py / 03_gene_map.py) ─────────
def read_gtf(path):
    """Read the gene and transcript records of a gzipped GTF file, in file order.

    Returns a gene table and a transcript table. Identifiers keep their version, since the
    chromosome Y copy of a pseudoautosomal gene differs from its chromosome X copy only in
    the suffix _PAR_Y after the version.
    """
    gene_rows = []
    transcript_rows = []
    with gzip.open(path, "rt") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, _, record_type, start, end, _, strand, _, attributes = fields
            if record_type not in ("gene", "transcript"):
                continue
            pairs = re.findall(r'(\w+) "([^"]*)"', attributes)
            values = dict(pairs)
            tags = [value for key, value in pairs if key == "tag"]
            if record_type == "gene":
                gene_rows.append(
                    {
                        "gene_id": values["gene_id"],
                        "symbol": values["gene_name"],
                        "hgnc_id": values.get("hgnc_id", ""),
                        "biotype": values["gene_type"],
                        "chrom": chrom,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                        "has_par_tag": "PAR" in tags,
                    }
                )
            else:
                transcript_rows.append(
                    {
                        "gene_id": values["gene_id"],
                        "transcript_type": values["transcript_type"],
                        "is_readthrough": "readthrough_transcript" in tags,
                    }
                )
    return pd.DataFrame(gene_rows), pd.DataFrame(transcript_rows)


def read_exons(path, gene_ids):
    """Read the exon records of the given (versioned) genes from a gzipped GTF file.

    gene_ids holds identifiers with their version, so the chromosome Y copy of a
    pseudoautosomal gene (identifier ends in _PAR_Y) is not read. Returns one row per exon
    of each transcript, with the gene's identifier, start and end.
    """
    rows = []
    with gzip.open(path, "rt") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] != "exon":
                continue
            gene_id = re.search(r'gene_id "([^"]*)"', fields[8]).group(1)
            if gene_id in gene_ids:
                rows.append({"gene_id": gene_id, "start": int(fields[3]), "end": int(fields[4])})
    return pd.DataFrame(rows)


def exon_overlap(exons_source, exons_target):
    """Fraction of the source gene's exon bases that lie inside the candidate's exons.

    Each table holds the exons of one gene, one row per exon of each transcript. A base is
    marked once however many exons cover it, so each gene counts the union of its exons.
    """
    first = min(exons_source["start"].min(), exons_target["start"].min())
    last = max(exons_source["end"].max(), exons_target["end"].max())
    in_source = np.zeros(last - first + 1, dtype=bool)
    in_target = np.zeros(last - first + 1, dtype=bool)
    for start, end in zip(exons_source["start"], exons_source["end"], strict=True):
        in_source[start - first : end - first + 1] = True
    for start, end in zip(exons_target["start"], exons_target["end"], strict=True):
        in_target[start - first : end - first + 1] = True
    return (in_source & in_target).sum() / in_source.sum()


def _genes_in(cell):
    """Return the semicolon-separated genes of a list cell, or [] for an empty cell."""
    return cell.split(";") if cell else []


def _gene_table(gtf_path, drop_par):
    """Read a GTF into a version-stripped gene table (gene_id, gene_version, …).

    drop_par removes the chromosome Y copies of pseudoautosomal genes before the version is
    stripped, so that the remaining identifiers are unique once versionless — used for the
    source table, whose identifiers are the (non-PAR) count-matrix columns.
    """
    genes, _ = read_gtf(gtf_path)
    if drop_par:
        genes = genes[~genes["has_par_tag"]].reset_index(drop=True)
    genes[["gene_id", "gene_version"]] = genes["gene_id"].str.split(".", expand=True)
    # Empty hgnc_id means "no HGNC identifier"; two missing identifiers never match.
    genes["hgnc_id"] = genes["hgnc_id"].replace("", np.nan)
    return genes


# ── Crosswalk computation (steps 1-6; pure, AnnData-independent) ──────────────────
def compute_gene_crosswalk(
    source_ids,
    source_gtf,
    target_gtf,
    stable_id_history,
    *,
    min_reciprocal_overlap=DEFAULT_MIN_RECIPROCAL_OVERLAP,
    min_exon_overlap=DEFAULT_MIN_EXON_OVERLAP,
    clone_based_name=CLONE_BASED_NAME,
    accepted_biotype_changes=None,
    checksums=None,
):
    """Decide, for each source gene, whether its column is kept, renamed or dropped.

    Parameters
    ----------
    source_ids : sequence of str
        Unversioned source Ensembl IDs in count-matrix order (``adata.var_names``).
    source_gtf, target_gtf, stable_id_history : path-like
        The source GENCODE GTF, the target GENCODE GTF, and the Ensembl stable-id history.
    accepted_biotype_changes : dict, optional
        ``{source_id: target_id}`` biotype changes that ``same_biotype`` accepts
        (the atlas supplies ``{ENSG00000184258: ENSG00000288642}`` for CDR1).
    checksums : dict, optional
        ``{"source_gtf": {"md5": …}, "target_gtf": {"md5": …},
        "stable_id_history": {"bsd_sum": …}}`` — verified before anything is read.

    Returns
    -------
    pandas.DataFrame
        One row per source gene, in input order, indexed 0..n-1, with at least:
        ``gene_id_source, symbol_source, biotype_source, action, gene_id_target,
        reason, in_object`` plus the per-retired-gene provenance (candidate, overlaps,
        the eight conditions).
    """
    accepted_biotype_changes = accepted_biotype_changes or {}
    source_ids = list(map(str, source_ids))

    # Step 1: verify checksums (no network).
    if checksums:
        verify_file(source_gtf, md5=checksums.get("source_gtf", {}).get("md5"))
        verify_file(target_gtf, md5=checksums.get("target_gtf", {}).get("md5"))
        verify_file(
            stable_id_history, bsd_sum=checksums.get("stable_id_history", {}).get("bsd_sum")
        )

    # Step 2: source gene table, restricted to the matrix columns, in matrix order.
    genes_src_all = _gene_table(source_gtf, drop_par=True).set_index("gene_id")
    missing = [gid for gid in source_ids if gid not in genes_src_all.index]
    if missing:
        sys.exit(
            f"FAIL {len(missing)} var_names absent from the source GTF "
            f"(e.g. {missing[:5]}): the source_gtf must match the counts' reference build"
        )
    genes_source = genes_src_all.loc[source_ids].reset_index()  # column 'gene_id', in order

    # Step 3: target gene table, version-stripped.
    genes_target = _gene_table(target_gtf, drop_par=False)
    genes_target["start"] = genes_target["start"].astype(int)
    genes_target["end"] = genes_target["end"].astype(int)
    target_by_id = genes_target.set_index("gene_id")

    source_id_set = set(source_ids)
    target_id_set = set(genes_target["gene_id"])
    kept_ids = source_id_set & target_id_set

    # Step 4: retired identifiers and their successors (02_retired_genes.py).
    retired = genes_source[~genes_source["gene_id"].isin(target_id_set)].reset_index(drop=True)

    #   First source: target genes with the same HGNC identifier.
    target_with_hgnc = genes_target[genes_target["hgnc_id"].notna()]
    target_ids_by_hgnc = target_with_hgnc.groupby("hgnc_id")["gene_id"].agg(";".join)
    retired["same_hgnc_id"] = retired["hgnc_id"].map(target_ids_by_hgnc)

    #   Second source: the Ensembl stable-id history.
    events = pd.read_csv(
        stable_id_history, sep="\t", header=None, names=EVENT_COLUMNS, dtype=str
    )
    successor_events = events[
        (events["type"] == "gene")
        & events["old_stable_id"].isin(retired["gene_id"])
        & events["new_stable_id"].isin(target_id_set)
    ]
    successor_pairs = successor_events[["old_stable_id", "new_stable_id"]].drop_duplicates()
    successor_pairs = successor_pairs.sort_values(["old_stable_id", "new_stable_id"])
    successors_by_id = successor_pairs.groupby("old_stable_id")["new_stable_id"].agg(";".join)
    retired["history_successors"] = retired["gene_id"].map(successors_by_id)

    #   Checks (02 lines 419-434).
    listed_genes = pd.concat(
        [
            retired["same_hgnc_id"].dropna().str.split(";").explode(),
            retired["history_successors"].dropna().str.split(";").explode(),
        ]
    )
    n_outside = int((~listed_genes.isin(target_id_set)).sum())
    if n_outside:
        sys.exit(f"FAIL {n_outside} listed successor genes outside the target build")
    n_several_hgnc = int(retired["same_hgnc_id"].str.contains(";", na=False).sum())
    if n_several_hgnc:
        sys.exit(
            f"FAIL {n_several_hgnc} retired identifiers with several target genes of their "
            "HGNC identifier (the candidate must be the one gene with the same HGNC id)"
        )

    # Step 5: choose the candidate, measure overlaps, test the eight conditions (03_gene_map.py).
    gene_map = pd.DataFrame(
        {
            "gene_id_source": retired["gene_id"],
            "symbol_source": retired["symbol"],
            "hgnc_id_source": retired["hgnc_id"],
            "biotype_source": retired["biotype"],
            "chrom_source": retired["chrom"],
            "start_source": retired["start"].astype(int),
            "end_source": retired["end"].astype(int),
            "strand_source": retired["strand"],
        }
    )
    same_hgnc_ids = retired["same_hgnc_id"]
    history_ids = retired["history_successors"]
    same_hgnc_genes = same_hgnc_ids.fillna("").map(_genes_in)
    history_genes = history_ids.fillna("").map(_genes_in)
    has_same_hgnc_gene = same_hgnc_ids.notna()
    has_single_history_gene = history_genes.map(len) == 1
    from_history = ~has_same_hgnc_gene & has_single_history_gene
    gene_map["candidate_target"] = same_hgnc_ids
    gene_map.loc[from_history, "candidate_target"] = history_ids[from_history]
    gene_map["candidate_source_kind"] = None
    gene_map.loc[has_same_hgnc_gene, "candidate_source_kind"] = "same HGNC identifier"
    gene_map.loc[from_history, "candidate_source_kind"] = "Ensembl history"
    for column in ["symbol", "hgnc_id", "biotype", "chrom", "start", "end", "strand"]:
        gene_map[f"{column}_target"] = gene_map["candidate_target"].map(target_by_id[column])

    #   The fraction of each span that the other covers.
    shared_bases = (
        np.minimum(gene_map["end_source"], gene_map["end_target"])
        - np.maximum(gene_map["start_source"], gene_map["start_target"])
        + 1
    ).clip(lower=0)
    elsewhere = gene_map["candidate_target"].notna() & (
        (gene_map["chrom_target"] != gene_map["chrom_source"])
        | (gene_map["strand_target"] != gene_map["strand_source"])
    )
    shared_bases[elsewhere] = 0
    gene_map["overlap_of_source"] = shared_bases / (
        gene_map["end_source"] - gene_map["start_source"] + 1
    )
    gene_map["overlap_of_candidate"] = shared_bases / (
        gene_map["end_target"] - gene_map["start_target"] + 1
    )

    #   The fraction of the source gene's exon bases inside the candidate's exons.
    versioned_source = retired["gene_id"] + "." + retired["gene_version"]
    versioned_target = (
        gene_map["candidate_target"]
        + "."
        + gene_map["candidate_target"].map(target_by_id["gene_version"])
    )
    exons_source = read_exons(source_gtf, set(versioned_source))
    exons_target = read_exons(target_gtf, set(versioned_target.dropna()))
    n_without_exons = len(set(versioned_source) - set(exons_source["gene_id"]))
    n_without_exons += len(set(versioned_target.dropna()) - set(exons_target["gene_id"]))
    if n_without_exons:
        sys.exit(f"FAIL {n_without_exons} retired genes or candidates without exon records")
    gene_map["exon_overlap_of_source"] = np.where(
        gene_map["candidate_target"].notna(), 0.0, np.nan
    )
    for row in gene_map.index[shared_bases > 0]:
        old_exons = exons_source[exons_source["gene_id"] == versioned_source[row]]
        new_exons = exons_target[exons_target["gene_id"] == versioned_target[row]]
        gene_map.loc[row, "exon_overlap_of_source"] = exon_overlap(old_exons, new_exons)

    #   The eight conditions.
    named_genes = same_hgnc_genes + history_genes
    claims = gene_map["candidate_target"].value_counts()
    gene_map["has_successor"] = named_genes.map(len) > 0
    gene_map["not_already_a_column"] = named_genes.map(kept_ids.isdisjoint)
    gene_map["single_candidate"] = has_same_hgnc_gene | has_single_history_gene
    gene_map["not_shared"] = gene_map["candidate_target"].map(claims) == 1
    same_symbol = gene_map["symbol_target"] == gene_map["symbol_source"]
    same_hgnc_id = gene_map["hgnc_id_target"] == gene_map["hgnc_id_source"]
    clone_based = (
        gene_map["candidate_target"].notna()
        & gene_map["hgnc_id_source"].isna()
        & gene_map["symbol_source"].str.fullmatch(clone_based_name)
    )
    gene_map["same_name"] = same_symbol | same_hgnc_id | clone_based
    same_biotype = gene_map["biotype_target"] == gene_map["biotype_source"]
    accepted_candidate = gene_map["gene_id_source"].map(accepted_biotype_changes)
    gene_map["same_biotype"] = same_biotype | (accepted_candidate == gene_map["candidate_target"])
    gene_map["same_locus"] = (gene_map["overlap_of_source"] >= min_reciprocal_overlap) & (
        gene_map["overlap_of_candidate"] >= min_reciprocal_overlap
    )
    gene_map["exons_in_candidate"] = gene_map["exon_overlap_of_source"] >= min_exon_overlap

    #   Decide: rename when every condition holds, else drop with the first failing reason.
    gene_map["reason"] = None
    for condition, reason in REASONS.items():
        first_failure = ~gene_map[condition] & gene_map["reason"].isna()
        gene_map.loc[first_failure, "reason"] = reason
    gene_map["action"] = np.where(gene_map["reason"].isna(), "rename", "drop")
    gene_map["gene_id_target"] = gene_map["candidate_target"].where(gene_map["action"] == "rename")

    #   Checks (03 lines 742-754).
    renamed = gene_map[gene_map["action"] == "rename"]
    if int(renamed["gene_id_target"].duplicated().sum()):
        sys.exit("FAIL rename targets named twice (two columns are never summed)")
    if int(renamed["gene_id_target"].isin(kept_ids).sum()):
        sys.exit("FAIL rename targets that are kept source identifiers (never summed)")

    # Step 6: full-width map over the source genes, in matrix order (08_var.py apply half).
    crosswalk = pd.DataFrame(
        {
            "gene_id_source": genes_source["gene_id"].to_numpy(),
            "symbol_source": genes_source["symbol"].to_numpy(),
            "biotype_source": genes_source["biotype"].to_numpy(),
        }
    )
    retired_detail = gene_map.set_index("gene_id_source")
    is_retired = crosswalk["gene_id_source"].isin(retired_detail.index)
    crosswalk["action"] = "keep"
    crosswalk["gene_id_target"] = crosswalk["gene_id_source"]
    crosswalk["reason"] = None
    crosswalk.loc[is_retired, "action"] = crosswalk.loc[is_retired, "gene_id_source"].map(
        retired_detail["action"]
    )
    crosswalk.loc[is_retired, "gene_id_target"] = crosswalk.loc[
        is_retired, "gene_id_source"
    ].map(retired_detail["gene_id_target"])
    crosswalk.loc[is_retired, "reason"] = crosswalk.loc[is_retired, "gene_id_source"].map(
        retired_detail["reason"]
    )
    crosswalk["in_object"] = crosswalk["action"] != "drop"

    #   Attach the per-retired-gene provenance (NaN for kept genes).
    provenance_columns = [
        "candidate_target",
        "candidate_source_kind",
        "symbol_target",
        "hgnc_id_target",
        "biotype_target",
        "overlap_of_source",
        "overlap_of_candidate",
        "exon_overlap_of_source",
        *REASONS.keys(),
    ]
    crosswalk = crosswalk.merge(
        retired_detail[provenance_columns],
        how="left",
        left_on="gene_id_source",
        right_index=True,
    )

    #   Uniqueness check over the kept columns (08 lines 861-867).
    in_object = crosswalk[crosswalk["in_object"]]
    if int(in_object["gene_id_target"].duplicated().sum()):
        sys.exit("FAIL repeated target identifiers in the object (the var index must be unique)")
    n_renamed_expected = int((gene_map["action"] == "rename").sum())
    n_relabelled = int((in_object["gene_id_source"] != in_object["gene_id_target"]).sum())
    if n_relabelled != n_renamed_expected:
        sys.exit(
            f"FAIL {n_relabelled} relabelled genes against {n_renamed_expected} renames in the map"
        )
    return crosswalk


# ── snakemake entry point (step 7: apply the crosswalk to the object) ─────────────
def _apply(adata, crosswalk, logger):
    """Relabel renamed columns, drop dropped columns, carry source provenance into var.

    Counts (X, layers, varm rows) are carried unchanged; two columns are never summed. In
    late mode (input = processed.h5ad) anndata's var-axis slicing also slices varm['PCs'].
    """
    source_ids = np.asarray(crosswalk["gene_id_source"].to_numpy(), dtype=object)
    if not np.array_equal(np.asarray(adata.var_names, dtype=object), source_ids):
        sys.exit("FAIL crosswalk gene_id_source does not match adata.var_names in order")

    # Carry provenance on the full-width var before subsetting.
    if "gene_name" in adata.var.columns:
        adata.var["gene_symbol_source"] = adata.var["gene_name"].to_numpy()
    else:
        adata.var["gene_symbol_source"] = crosswalk["symbol_source"].to_numpy()
    adata.var["gene_id_source"] = source_ids
    adata.var["feature_is_filtered"] = False

    keep = crosswalk["in_object"].to_numpy()
    target = crosswalk["gene_id_target"].to_numpy()
    adata = adata[:, keep].copy()
    adata.var_names = pd.Index(target[keep], name=None)

    if not adata.var_names.is_unique:
        sys.exit("FAIL reconciled var_names are not unique")
    n_keep = int((crosswalk["action"] == "keep").sum())
    n_rename = int((crosswalk["action"] == "rename").sum())
    n_drop = int((crosswalk["action"] == "drop").sum())
    logger.info(f"Gene lift: {n_keep} kept, {n_rename} renamed, {n_drop} dropped")
    logger.info(f"Object now carries {adata.n_vars} genes (in_object) of {len(crosswalk)} source")
    return adata


def _main():
    import logging

    import anndata as ad

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
    )
    logger = logging.getLogger(__name__)

    def _excepthook(exc_type, exc_value, exc_tb):
        logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_tb))
        sys.__excepthook__(exc_type, exc_value, exc_tb)

    sys.excepthook = _excepthook

    input_h5ad = snakemake.input.object
    source_gtf = snakemake.input.source_gtf
    target_gtf = snakemake.input.target_gtf
    stable_id_history = snakemake.input.stable_id_history
    output_h5ad = snakemake.output.reconciled
    output_crosswalk = snakemake.output.crosswalk

    params = dict(snakemake.params.get("reconcile", {}) or {})
    accepted = params.get("accepted_biotype_changes", {}) or {}
    checksums = {
        key: {ck: ref[ck] for ck in ("md5", "bsd_sum") if ck in ref}
        for key, ref in (params.get("references", {}) or {}).items()
    }

    logger.info(f"Reconciling gene space of: {input_h5ad}")
    adata = ad.read_h5ad(input_h5ad)
    logger.info(f"Input shape: {adata.n_obs} cells x {adata.n_vars} genes")

    crosswalk = compute_gene_crosswalk(
        list(adata.var_names),
        source_gtf,
        target_gtf,
        stable_id_history,
        min_reciprocal_overlap=params.get("min_reciprocal_overlap", DEFAULT_MIN_RECIPROCAL_OVERLAP),
        min_exon_overlap=params.get("min_exon_overlap", DEFAULT_MIN_EXON_OVERLAP),
        accepted_biotype_changes=accepted,
        checksums=checksums or None,
    )

    Path(output_crosswalk).parent.mkdir(parents=True, exist_ok=True)
    crosswalk.to_csv(output_crosswalk, index=False)
    logger.info(f"Wrote gene crosswalk: {output_crosswalk}")

    adata = _apply(adata, crosswalk, logger)

    logger.info(f"Saving reconciled object to: {output_h5ad}")
    Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad)
    logger.info("Done.")


if "snakemake" in globals():
    _main()
