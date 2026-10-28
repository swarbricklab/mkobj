#!/usr/bin/env python3
"""
workflow/scripts/resolve_keep.py
Stage 2 of the config-gated keep-list filter (see rules/filter_deposit.smk).

Declaratively resolve which cells the deposit keeps, with DuckDB. The obs CSV
(stage 1) is LEFT JOINed against a *disposition* table — a projection of the
project samplesheet that maps a resolution key (sample_id / tissue_id) to a
usage/disposition value — and a config predicate over (obs columns x the joined
value column) decides the keep-list.

Everything is config-driven (nothing hard-coded to one deposit):

    filter:
      disposition_table: path to the disposition CSV (sample_id,tissue_id,usage,...)
      disposition_filter: {assay: chromium}        # optional row pre-filter (AND of equalities)
      disposition_value_col: usage                 # the column joined onto obs
      joins:                                        # (obs_col, disposition_col) probed in order
        - [sample_id, sample_id]
        - [sample_id, tissue_id]
      aliases: {"4270": "4270-1", ...}             # applied to obs join-column values first
      exclude_capture_regex: buffy                  # optional; dropped regardless of predicate
      exclude_column: capture_id
      predicate: "status = 'singlet' AND usage IN ('analysed','control')"

Outputs (all tracked, so they document the deposit's inclusion criteria):
  keep_cells.csv       one column `cell_id` — the cells the subset keeps
  dropped_cells.csv    cell_id + status + value + exclude col + `reason`
  status_usage_breakdown.csv   status x value x kept, with cell counts

Design notes / gotchas handled:
  * Join key is the UNIQUE capture-prefixed cell_id from obs (not the 16bp barcode).
  * LEFT JOIN: unmatched rows get NULL value and are dropped with reason
    'unmatched_usage' — the NULL is handled and documented, never silently kept.
  * The disposition projection is deduped to 1 row per key; keys that map to >1
    distinct value are reported (conflicts) and collapsed deterministically.
  * `joins` is an ordered COALESCE: the first disposition column that matches
    the (aliased) obs value wins — e.g. sample_id first, tissue_id as fallback.
"""

import sys
import logging

import duckdb

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)


def q(ident: str) -> str:
    """Quote a SQL identifier."""
    return '"' + str(ident).replace('"', '""') + '"'


def lit(value) -> str:
    """Quote a SQL string literal."""
    return "'" + str(value).replace("'", "''") + "'"


# ── Inputs / params ──────────────────────────────────────────────────────────────
obs_csv = snakemake.input.obs_csv
disposition_csv = snakemake.input.disposition

keep_out = snakemake.output.keep
dropped_out = snakemake.output.dropped
breakdown_out = snakemake.output.breakdown

value_col = snakemake.params.get("disposition_value_col", "usage")
joins = snakemake.params.get("joins") or [["sample_id", "sample_id"]]
aliases = snakemake.params.get("aliases") or {}
disposition_filter = snakemake.params.get("disposition_filter") or {}
predicate = snakemake.params["predicate"]
exclude_regex = snakemake.params.get("exclude_capture_regex") or None
exclude_col = snakemake.params.get("exclude_column", "capture_id")

logger.info(f"obs CSV: {obs_csv}")
logger.info(f"disposition CSV: {disposition_csv}")
logger.info(f"joins (obs_col -> disposition_col): {joins}")
logger.info(f"value column: {value_col}; predicate: {predicate}")
logger.info(f"aliases: {aliases}; exclude {exclude_col} =~ {exclude_regex!r}")

con = duckdb.connect()

# Read both tables as all-varchar so types never surprise the join/predicate;
# empty strings become NULL.
con.execute(
    "CREATE TABLE obs AS "
    "SELECT * FROM read_csv_auto(?, all_varchar=true, nullstr='')", [obs_csv]
)
con.execute(
    "CREATE TABLE disp_raw AS "
    "SELECT * FROM read_csv_auto(?, all_varchar=true, nullstr='')", [disposition_csv]
)
n_obs = con.execute("SELECT count(*) FROM obs").fetchone()[0]
logger.info(f"Loaded obs rows: {n_obs}")

# ── Disposition pre-filter (optional AND of equalities) ──────────────────────────
where = " AND ".join(f"{q(c)} = {lit(v)}" for c, v in disposition_filter.items())
disp_where = f"WHERE {where}" if where else ""
con.execute(f"CREATE TABLE disp AS SELECT * FROM disp_raw {disp_where}")
n_disp = con.execute("SELECT count(*) FROM disp").fetchone()[0]
logger.info(f"Disposition rows after filter ({disposition_filter or 'none'}): {n_disp}")

# ── Aliased obs join columns ─────────────────────────────────────────────────────
obs_join_cols = sorted({j[0] for j in joins})
alias_exprs = []
for col in obs_join_cols:
    if aliases:
        cases = " ".join(f"WHEN {lit(k)} THEN {lit(v)}" for k, v in aliases.items())
        expr = f"CASE {q(col)} {cases} ELSE {q(col)} END"
    else:
        expr = q(col)
    alias_exprs.append(f"{expr} AS {q('__ak_' + col)}")
con.execute(f"CREATE TABLE obs_a AS SELECT *, {', '.join(alias_exprs)} FROM obs")

# ── One deduped key->value map per disposition column, + conflict report ─────────
disp_cols = []
for _, dcol in joins:
    if dcol not in disp_cols:
        disp_cols.append(dcol)
for i, dcol in enumerate(disp_cols):
    con.execute(
        f"CREATE TABLE map_{i} AS "
        f"SELECT {q(dcol)} AS k, "
        f"       any_value({q(value_col)}) AS v, "
        f"       count(DISTINCT {q(value_col)}) AS n_distinct "
        f"FROM disp WHERE {q(dcol)} IS NOT NULL GROUP BY {q(dcol)}"
    )
    conflicts = con.execute(
        f"SELECT k, n_distinct FROM map_{i} WHERE n_distinct > 1"
    ).fetchall()
    if conflicts:
        logger.warning(
            f"disposition column {dcol!r}: {len(conflicts)} key(s) map to >1 "
            f"distinct {value_col} — collapsed deterministically: {conflicts[:10]}"
        )

# ── Resolve value via ordered COALESCE over the joins ────────────────────────────
coalesce_terms = []
join_clauses = []
for n, (ocol, dcol) in enumerate(joins):
    i = disp_cols.index(dcol)
    coalesce_terms.append(f"m{n}.v")
    join_clauses.append(
        f"LEFT JOIN map_{i} m{n} ON obs_a.{q('__ak_' + ocol)} = m{n}.k"
    )
value_expr = f"COALESCE({', '.join(coalesce_terms)})"

con.execute(
    f"CREATE TABLE resolved AS "
    f"SELECT obs_a.*, {value_expr} AS {q(value_col)} "
    f"FROM obs_a {' '.join(join_clauses)}"
)

# ── Keep / drop classification ───────────────────────────────────────────────────
if exclude_regex:
    excl = (
        f"regexp_matches(lower(CAST({q(exclude_col)} AS VARCHAR)), "
        f"lower({lit(exclude_regex)}))"
    )
else:
    excl = "FALSE"

# keep = NOT excluded AND predicate true. COALESCE to FALSE so an unmatched row
# (value IS NULL -> predicate evaluates to SQL NULL) is a definite drop, never a
# three-valued NULL that slips past both the keep and the drop filters. reason
# documents every drop, with the LEFT-JOIN NULL (unmatched) called out before a
# generic predicate failure.
keep_flag = f"COALESCE((NOT ({excl})) AND ({predicate}), FALSE)"
reason = (
    f"CASE "
    f"WHEN {excl} THEN 'excluded_capture' "
    f"WHEN {q(value_col)} IS NULL THEN 'unmatched_usage' "
    f"WHEN NOT ({predicate}) THEN 'predicate_fail' "
    f"ELSE NULL END"
)
con.execute(
    f"CREATE TABLE classified AS "
    f"SELECT *, ({keep_flag}) AS keep_flag, ({reason}) AS reason FROM resolved"
)

n_keep = con.execute("SELECT count(*) FROM classified WHERE keep_flag").fetchone()[0]
n_drop = n_obs - n_keep
logger.info(f"KEEP {n_keep} / DROP {n_drop} (of {n_obs})")
for r, c in con.execute(
    "SELECT reason, count(*) FROM classified WHERE NOT keep_flag "
    "GROUP BY reason ORDER BY 2 DESC"
).fetchall():
    logger.info(f"  drop reason {r}: {c}")

# ── Emit artifacts ───────────────────────────────────────────────────────────────
con.execute(
    f"COPY (SELECT cell_id FROM classified WHERE keep_flag ORDER BY cell_id) "
    f"TO {lit(keep_out)} (HEADER, DELIMITER ',')"
)
con.execute(
    f"COPY (SELECT cell_id, status, {q(value_col)}, {q(exclude_col)}, reason "
    f"FROM classified WHERE NOT keep_flag ORDER BY reason, cell_id) "
    f"TO {lit(dropped_out)} (HEADER, DELIMITER ',')"
)
con.execute(
    f"COPY (SELECT status, COALESCE({q(value_col)}, 'UNRESOLVED') AS {q(value_col)}, "
    f"keep_flag AS kept, count(*) AS n_cells FROM classified "
    f"GROUP BY status, {q(value_col)}, keep_flag ORDER BY n_cells DESC) "
    f"TO {lit(breakdown_out)} (HEADER, DELIMITER ',')"
)
logger.info(f"Wrote keep={keep_out} dropped={dropped_out} breakdown={breakdown_out}")
con.close()
