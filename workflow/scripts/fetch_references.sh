#!/usr/bin/env bash
# Fetch the three gene references the reconcile_gene_space stage needs, at pinned versions, and
# check each against its published checksum. Generalised from the CELLxGENE curation step
# 00_download_references.sh of the publication repo.
#
# Environment: none. Uses curl, md5sum and sum, and needs the network. NCI login/ps nodes and the
# copyq data-mover have it; compute nodes do NOT. Run this on a data-mover, e.g.
#   qx --internet --env dt3 --storage gdata/a56+scratch/a56 -P a56 -- \
#     bash modules/mkobj/workflow/scripts/fetch_references.sh resources/gene_reference
#
# The three pins below mirror config/mkobj/config.yaml -> reconcile.references. The compute stage
# (reconcile_gene_space.py, step 1) re-verifies the same checksums with no network.
#
# Argument 1 (optional): output directory. Default: resources/gene_reference (repo-root relative).
#
# Outputs, in the output directory:
# - gencode.v32.primary_assembly.annotation.gtf.gz        (GENCODE release 32; 2020-A gene set)
# - gencode.v48.primary_assembly.annotation.gtf.gz        (GENCODE release 48; CxG 7.1.0 gene set)
# - homo_sapiens_core_114_38_stable_id_event.txt.gz       (Ensembl 114 stable-id history)
#
# A checksum that differs stops the script.
set -euo pipefail

directory="${1:-resources/gene_reference}"
gencode=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human
ensembl=https://ftp.ensembl.org/pub/release-114/mysql/homo_sapiens_core_114_38

# Print the check, and stop when the observed value differs from the expected one.
check() {  # check <quantity> <observed> <expected>
  if [[ "$2" != "$3" ]]; then
    echo "FAIL $1: $2 against $3 (published checksum)" >&2
    exit 1
  fi
  echo "ok   $1: $2"
}

mkdir -p "$directory"

# 1. GENCODE v32, the gene set of the Cell Ranger reference refdata-gex-GRCh38-2020-A
file=gencode.v32.primary_assembly.annotation.gtf.gz
curl -fsSL -o "$directory/$file" "$gencode/release_32/$file"
check "md5 of $file" "$(md5sum < "$directory/$file" | cut -d " " -f 1)" \
  f660ce967d954fe704700ddc81ebf373

# 2. GENCODE v48, the gene set CELLxGENE schema 7.1.0 requires
file=gencode.v48.primary_assembly.annotation.gtf.gz
curl -fsSL -o "$directory/$file" "$gencode/release_48/$file"
check "md5 of $file" "$(md5sum < "$directory/$file" | cut -d " " -f 1)" \
  1b68ecb47fd4cc3946819fec58a59d8f

# 3. Ensembl 114 stable identifier history
file=homo_sapiens_core_114_38_stable_id_event.txt.gz
curl -fsSL -o "$directory/$file" "$ensembl/stable_id_event.txt.gz"
check "BSD checksum of $file" "$(sum < "$directory/$file")" "22726 16128"

# 4. CELLxGENE's own accepted-gene list (genes_homo_sapiens.csv.gz) — the exact set
# `cellxgene-schema validate` checks. Copied out of the cxg_schema conda env so it always
# matches the pinned validator version (no GitHub tag/md5 to guess). Best-effort: if the
# env/package is unavailable, reconcile_gene_space.py skips the early membership check and
# validate_cellxgene backstops it, so this does not fail the fetch.
file=genes_homo_sapiens.csv.gz
source /g/data/a56/conda/miniconda/etc/profile.d/conda.sh 2>/dev/null || true
if src=$(conda run -n cxg_schema python -c \
      "import cellxgene_schema,os;print(os.path.join(os.path.dirname(cellxgene_schema.__file__),'gencode_files','$file'))" \
      2>/dev/null) && [[ -f "$src" ]]; then
  cp "$src" "$directory/$file"
  echo "ok   $file copied from cxg_schema env ($(gzip -dc "$directory/$file" | grep -c '^ENSG') gene ids)"
else
  echo "warn $file NOT sourced (cxg_schema env/cellxgene_schema unavailable); the reconcile" \
       "membership check will be skipped — validate_cellxgene backstops it" >&2
fi

echo "done: references in $directory"
