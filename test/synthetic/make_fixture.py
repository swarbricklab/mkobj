#!/usr/bin/env python3
"""
test/synthetic/make_fixture.py
Generate a small, fully synthetic mkobj input set that runs on a laptop.

The DVC-tracked test data (config/test.yaml) lives on the cluster. This fixture
exercises the same code paths locally, with known ground truth:

  * 2 captures x 2 donors each (S1, S2 in CAP_A; S3, S4 in CAP_B).
  * 4 fine cell types in 2 major lineages (Epithelial: LumHR, Basal;
    T-cells: CD4, CD8) with distinct marker programmes.
  * ~8% doublets made by summing two cells' counts. Cross-donor doublets are
    called 'doublet' by the fake demux; within-donor ones look like singlets.
  * samples.csv lists S1, S2, S3 only: S4 is a passenger sample, so cohort
    subsetting has something to drop.

Outputs (under --out, default test/synthetic/data):
  cellranger/<capture>/{barcodes,features}.tsv.gz, matrix.mtx.gz
  demux/<capture>/cell_assignment.tsv, sample_barcode_assignments.tsv
  annotation/<capture>/cell_types.csv      (cellranger-annotate-style labels)
  lineage_crosswalk.csv                    (coarse label -> major lineage)
  truth/<capture>.csv                      (barcode, donor(s), fine/major type, is_doublet)
  captures.csv, samples.csv

Usage:
  python test/synthetic/make_fixture.py [--out DIR] [--seed N]
"""

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, sparse

CAPTURES = {"CAP_A": ["S1", "S2"], "CAP_B": ["S3", "S4"]}
COHORT = ["S1", "S2", "S3"]
DONOR_SEX = {"S1": "F", "S2": "M", "S3": "F", "S4": "M"}

FINE_TO_MAJOR = {
    "LumHR": "Epithelial",
    "Basal": "Epithelial",
    "CD4": "T-cells",
    "CD8": "T-cells",
}
# Marker genes per fine type (shared major markers listed under both).
MARKERS = {
    "LumHR": ["EPCAM", "KRT8", "KRT18", "ESR1", "PGR", "AR", "FOXA1", "GATA3"],
    "Basal": ["EPCAM", "KRT5", "KRT14", "KRT17", "ACTA2", "MYLK", "TP63", "OXTR"],
    "CD4": ["PTPRC", "CD3E", "CD3D", "CD2", "CD4", "IL7R", "CCR7", "FOXP3"],
    "CD8": ["PTPRC", "CD3E", "CD3D", "CD2", "CD8A", "CD8B", "GZMK", "NKG7"],
}
# cellranger annotate-style labels (coarse_cell_type, fine_cell_type).
CR_LABEL = {
    "LumHR": ("epithelial cell", "luminal epithelial cell of mammary gland"),
    "Basal": ("epithelial cell", "basal cell"),
    "CD4": ("T cell", "CD4-positive, alpha-beta T cell"),
    "CD8": ("T cell", "CD8-positive, alpha-beta T cell"),
}
# Crosswalk from coarse label to the synthetic model's major vocabulary
# ("unknown cell type" is deliberately left unmapped).
CROSSWALK = {"epithelial cell": "Epithelial", "T cell": "T-cells"}
SEX_GENES = {"F": ["XIST"], "M": ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY"]}
HLA = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRA"]
N_FILLER = 900
N_CELLS = 700           # per capture, before doublets
DOUBLET_RATE = 0.08


def gene_table():
    named = sorted({g for gs in MARKERS.values() for g in gs})
    named += ["XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY"] + HLA
    named += [f"MT-CO{i}" for i in (1, 2, 3)] + ["RPS3", "RPL13", "RPS27"]
    names = named + [f"GENE{i:04d}" for i in range(N_FILLER)]
    ids = [f"ENSG{90000000000 + i:011d}" for i in range(len(names))]
    return pd.DataFrame({"gene_id": ids, "gene_name": names})


def simulate_cell(rng, genes, fine, donor, depth):
    """Poisson-ish counts from a cell-type + donor programme."""
    name_idx = {n: i for i, n in enumerate(genes.gene_name)}
    mu = rng.gamma(0.4, 1.0, size=len(genes))          # background
    for g in MARKERS[fine]:
        mu[name_idx[g]] += 25.0
    for g in SEX_GENES[DONOR_SEX[donor]]:
        mu[name_idx[g]] += 8.0
    for i, g in enumerate(HLA):                          # donor-specific HLA level
        mu[name_idx[g]] += 3.0 + 4.0 * ((int(donor[1:]) + i) % 3)
    for g in ("MT-CO1", "MT-CO2", "MT-CO3"):
        mu[name_idx[g]] += 6.0
    for g in ("RPS3", "RPL13", "RPS27"):
        mu[name_idx[g]] += 10.0
    mu = mu / mu.sum() * depth
    return rng.poisson(mu)


def make_capture(rng, capture, donors, genes, out):
    fine_types = list(FINE_TO_MAJOR)
    cells, truth = [], []
    for _ in range(N_CELLS):
        donor = donors[rng.integers(len(donors))]
        fine = fine_types[rng.integers(len(fine_types))]
        depth = rng.lognormal(np.log(2500), 0.35)
        cells.append(simulate_cell(rng, genes, fine, donor, depth))
        truth.append(dict(donor=donor, donor2="", fine=fine, fine2="", is_doublet=False))

    n_dbl = int(round(N_CELLS * DOUBLET_RATE))
    for _ in range(n_dbl):
        i, j = rng.choice(N_CELLS, size=2, replace=False)
        cells.append(cells[i] + cells[j])
        ti, tj = truth[i], truth[j]
        truth.append(dict(donor=ti["donor"], donor2=tj["donor"], fine=ti["fine"],
                          fine2=tj["fine"], is_doublet=True))

    order = rng.permutation(len(cells))
    X = np.vstack([cells[k] for k in order])
    truth = pd.DataFrame([truth[k] for k in order])
    barcodes = [f"{''.join(rng.choice(list('ACGT'), 16))}-1" for _ in range(len(order))]
    truth.insert(0, "barcode", barcodes)
    truth["major"] = truth["fine"].map(FINE_TO_MAJOR)

    # Cell Ranger filtered matrix (features x cells, gzip)
    cr = out / "cellranger" / capture
    cr.mkdir(parents=True, exist_ok=True)
    with gzip.open(cr / "matrix.mtx.gz", "wb") as fh:
        io.mmwrite(fh, sparse.csr_matrix(X.T).astype(np.int64), field="integer")
    feats = genes.assign(feature_type="Gene Expression")
    feats.to_csv(cr / "features.tsv.gz", sep="\t", header=False, index=False,
                 compression="gzip")
    pd.Series(barcodes).to_csv(cr / "barcodes.tsv.gz", sep="\t", header=False,
                               index=False, compression="gzip")

    # Fake SNP demux: cross-donor doublets -> doublet; 2% singlets unassigned.
    cross = truth["is_doublet"] & (truth["donor"] != truth["donor2"])
    unassigned = (~truth["is_doublet"]) & (rng.random(len(truth)) < 0.02)
    singlet = ~cross & ~unassigned
    dm = out / "demux" / capture
    dm.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "barcode": barcodes,
        "assignment": np.where(singlet, truth["donor"], "NA"),
        "status": np.where(singlet, "singlet", "NA"),
    }).to_csv(dm / "cell_assignment.tsv", sep="\t", index=False)
    pd.DataFrame({
        "Barcode": barcodes,
        "AtLeastHalfSinglet_DropletType": np.where(cross, "doublet", "singlet"),
    }).to_csv(dm / "sample_barcode_assignments.tsv", sep="\t", index=False)

    # cellranger annotate-style labels (with one unmapped oddball label)
    an = out / "annotation" / capture
    an.mkdir(parents=True, exist_ok=True)
    coarse = truth["fine"].map(lambda f: CR_LABEL[f][0])
    fine = truth["fine"].map(lambda f: CR_LABEL[f][1])
    odd = rng.random(len(coarse)) < 0.01
    coarse[odd] = fine[odd] = "unknown cell type"
    pd.DataFrame({"barcode": barcodes, "coarse_cell_type": coarse,
                  "fine_cell_type": fine}).to_csv(an / "cell_types.csv", index=False)

    (out / "truth").mkdir(parents=True, exist_ok=True)
    truth.assign(capture_id=capture, demux_doublet=cross.values).to_csv(
        out / "truth" / f"{capture}.csv", index=False)
    return len(barcodes)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--out", default=str(Path(__file__).parent / "data"))
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    out = Path(args.out)
    rng = np.random.default_rng(args.seed)
    genes = gene_table()

    counts = {c: make_capture(rng, c, d, genes, out) for c, d in CAPTURES.items()}
    pd.DataFrame({
        "capture": list(CAPTURES),
        "pool_id": ["POOL_1", "POOL_1"],
        "chip": ["standard", "ht"],
    }).to_csv(out / "captures.csv", index=False)
    rows = [dict(sample_id=s, capture_id=c, donor_id=f"D{s[1:]}", tissue_id=f"T{s[1:]}")
            for c, ss in CAPTURES.items() for s in ss if s in COHORT]
    pd.DataFrame(rows).to_csv(out / "samples.csv", index=False)
    pd.DataFrame({"cellranger_label": list(CROSSWALK), "lineage_major": list(CROSSWALK.values())}
                 ).to_csv(out / "lineage_crosswalk.csv", index=False)
    print(f"Wrote synthetic fixture to {out}: {counts}")


if __name__ == "__main__":
    main()
