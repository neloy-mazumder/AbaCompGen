#!/usr/bin/env bash
# =============================================================================
# STEP 02 — Multi-Locus Sequence Typing (MLST)
# conda env : clinical_env
# Tool      : mlst (Seemann) — wraps PubMLST databases
# Paper ref : "Pasteur scheme — 7 housekeeping genes"
# Runs BOTH Pasteur AND Oxford schemes for cross-validation.
# Output    : TSV with ST per genome for both schemes
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

OUT="$RESULTS_DIR/02_mlst"
LOG="$RESULTS_DIR/00_logs/02_mlst.log"
mkdir -p "$OUT"

echo "[INFO] Starting MLST — $(date)" | tee "$LOG"

# ─── PASTEUR SCHEME ──────────────────────────────────────────────────────────
echo "[INFO] Running Pasteur scheme (abaumannii_2)..." | tee -a "$LOG"

conda run -n clinical_env mlst \
    --scheme abaumannii_2 \
    --threads "$THREADS" \
    "$ASSEMBLY_DIR"/*.fna \
    > "$OUT/mlst_pasteur_raw.tsv" 2>>"$LOG"

# ─── OXFORD SCHEME ───────────────────────────────────────────────────────────
echo "[INFO] Running Oxford scheme (abaumannii)..." | tee -a "$LOG"

conda run -n clinical_env mlst \
    --scheme abaumannii \
    --threads "$THREADS" \
    "$ASSEMBLY_DIR"/*.fna \
    > "$OUT/mlst_oxford_raw.tsv" 2>>"$LOG"

# ─── FORMAT OUTPUT TABLES ────────────────────────────────────────────────────
echo "[INFO] Formatting output tables..." | tee -a "$LOG"

HEADER_PASTEUR="sample\tscheme\tST\tcpn60\tfusA\tgltA\tpyrG\trecA\trplB\trpoB"
HEADER_OXFORD="sample\tscheme\tST\tgltA\tgyrB\tgdhB\trecA\tcpn60\tgpi\trpoD"

for SCHEME in pasteur oxford; do
    RAW="$OUT/mlst_${SCHEME}_raw.tsv"
    CLEAN="$OUT/mlst_${SCHEME}.tsv"

    if [[ "$SCHEME" == "pasteur" ]]; then
        echo -e "$HEADER_PASTEUR" > "$CLEAN"
    else
        echo -e "$HEADER_OXFORD" > "$CLEAN"
    fi

    awk 'BEGIN{OFS="\t"} {
        n = split($1, a, "/"); name = a[n]
        sub(/\.fna$/, "", name)
        print name, $2, $3, $4, $5, $6, $7, $8, $9, $10
    }' "$RAW" >> "$CLEAN"

    echo "[INFO] Scheme $SCHEME formatted → $CLEAN" | tee -a "$LOG"
done

# ─── MERGE SCHEMES ───────────────────────────────────────────────────────────
MERGE="$RESULTS_DIR/13_summary_tables/mlst_combined.tsv"
echo -e "sample\tST_pasteur\tST_oxford" > "$MERGE"

join -t $'\t' -1 1 -2 1 \
    <(awk -F'\t' 'NR>1{print $1"\t"$3}' "$OUT/mlst_pasteur.tsv" | sort -k1,1) \
    <(awk -F'\t' 'NR>1{print $1"\t"$3}' "$OUT/mlst_oxford.tsv"  | sort -k1,1) \
    >> "$MERGE" 2>>"$LOG" || true

echo ""
echo "=== ST Frequency (Pasteur scheme) ===" | tee -a "$LOG"
awk -F'\t' 'NR>1 && $3 != "-" {print $3}' "$OUT/mlst_pasteur.tsv" \
    | sort | uniq -c | sort -rn \
    | awk '{printf "  ST%-8s  n=%d\n", $2, $1}' | tee -a "$LOG"

echo ""
echo "[DONE] Step 02 complete — $(date)" | tee -a "$LOG"
