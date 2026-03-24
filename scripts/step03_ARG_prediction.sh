#!/usr/bin/env bash
# =============================================================================
# STEP 03 — Antibiotic Resistance Gene (ARG) Prediction
# conda env : clinical_env
# Tools     : abricate (CARD db)  — mirrors paper's RGI approach
#             amrfinderplus (NCBI) — broader, well-curated complement
#             abricate (ResFinder) — acquired resistance genes
# Paper ref : "RGI tool / CARD database — perfect + strict criteria"
# Output    : Per-sample TSV + merged presence/absence matrix
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

OUT="$RESULTS_DIR/03_ARG_RGI"
OUT_AMR="$RESULTS_DIR/04_ARG_amrfinder"
LOG="$RESULTS_DIR/00_logs/03_ARG.log"
mkdir -p "$OUT" "$OUT_AMR"

echo "[INFO] Starting ARG prediction — $(date)" | tee "$LOG"

# ─── PART A: Abricate + CARD ──────────────────────────────────────────────────
echo "[INFO] === PART A: Abricate + CARD ===" | tee -a "$LOG"

conda run -n clinical_env abricate \
    --db card \
    --threads "$THREADS" \
    --minid 80 \
    --mincov 80 \
    "$ASSEMBLY_DIR"/*.fna \
    > "$OUT/abricate_card_raw.tsv" 2>>"$LOG"

awk -F'\t' 'NR>1 {
    n=split($1,a,"/"); name=a[n]; sub(/\.fna$/,"",name)
    print > "'"$OUT"'/" name "_card.tsv"
}' "$OUT/abricate_card_raw.tsv"

conda run -n clinical_env abricate \
    --summary "$OUT/abricate_card_raw.tsv" \
    > "$OUT/abricate_card_summary.tsv" 2>>"$LOG"

echo "[INFO] CARD summary → $OUT/abricate_card_summary.tsv" | tee -a "$LOG"

# ─── PART B: AMRFinderPlus ────────────────────────────────────────────────────
echo "[INFO] === PART B: AMRFinderPlus ===" | tee -a "$LOG"

AMR_ALL="$OUT_AMR/amrfinder_all_samples.tsv"
FIRST_SAMPLE=$(head -1 "$SAMPLE_LIST")

# Write header from first sample
conda run -n clinical_env amrfinder \
    --nucleotide "$ASSEMBLY_DIR/${FIRST_SAMPLE}.fna" \
    --organism Acinetobacter_baumannii \
    --plus \
    --threads "$THREADS" \
    --name "$FIRST_SAMPLE" \
    --database "$AMRFINDER_DB" \
    2>>"$LOG" | head -1 > "$AMR_ALL"

while IFS= read -r SAMPLE; do
    GENOME="$ASSEMBLY_DIR/${SAMPLE}.fna"
    [[ ! -f "$GENOME" ]] && continue

    SAMPLE_OUT="$OUT_AMR/${SAMPLE}_amrfinder.tsv"
    echo "[RUN ] AMRFinder: $SAMPLE" | tee -a "$LOG"

    conda run -n clinical_env amrfinder \
        --nucleotide "$GENOME" \
        --organism Acinetobacter_baumannii \
        --plus \
        --threads "$THREADS" \
        --name "$SAMPLE" \
        --output "$SAMPLE_OUT" \
        --database "$AMRFINDER_DB" \
        2>>"$LOG" || \
    conda run -n clinical_env amrfinder \
        --nucleotide "$GENOME" \
        --plus \
        --threads "$THREADS" \
        --name "$SAMPLE" \
        --output "$SAMPLE_OUT" \
        --database "$AMRFINDER_DB" \
        2>>"$LOG"

    tail -n +2 "$SAMPLE_OUT" >> "$AMR_ALL" 2>>"$LOG" || true
    echo "[DONE] $SAMPLE" | tee -a "$LOG"
done < "$SAMPLE_LIST"

# ─── PART C: Abricate + ResFinder ────────────────────────────────────────────
echo "[INFO] === PART C: Abricate + ResFinder DB ===" | tee -a "$LOG"

RESF_OUT="$RESULTS_DIR/14_abricate_resfinder"
mkdir -p "$RESF_OUT"

conda run -n clinical_env abricate \
    --db resfinder \
    --threads "$THREADS" \
    --minid 90 \
    --mincov 60 \
    "$ASSEMBLY_DIR"/*.fna \
    > "$RESF_OUT/abricate_resfinder_raw.tsv" 2>>"$LOG"

conda run -n clinical_env abricate \
    --summary "$RESF_OUT/abricate_resfinder_raw.tsv" \
    > "$RESF_OUT/abricate_resfinder_summary.tsv" 2>>"$LOG"

# ─── SUMMARY COUNT TABLE ─────────────────────────────────────────────────────
echo "[INFO] Generating ARG count table..." | tee -a "$LOG"
COUNT_TABLE="$RESULTS_DIR/13_summary_tables/ARG_counts.tsv"
echo -e "sample\tn_ARG_card\tn_ARG_amrfinder" > "$COUNT_TABLE"

while IFS= read -r SAMPLE; do
    N_CARD=$(grep -c "/${SAMPLE}\.fna" "$OUT/abricate_card_raw.tsv" 2>/dev/null || echo 0)
    N_AMR=0
    [[ -f "$OUT_AMR/${SAMPLE}_amrfinder.tsv" ]] && \
        N_AMR=$(awk 'NR>1' "$OUT_AMR/${SAMPLE}_amrfinder.tsv" | wc -l)
    echo -e "${SAMPLE}\t${N_CARD}\t${N_AMR}" >> "$COUNT_TABLE"
done < "$SAMPLE_LIST"

echo "[INFO] ARG count table → $COUNT_TABLE"
echo "[DONE] Step 03 complete — $(date)" | tee -a "$LOG"

echo ""
echo "=== OXA-type carbapenemase distribution ===" | tee -a "$LOG"
grep -i "OXA" "$OUT/abricate_card_raw.tsv" 2>/dev/null \
    | awk -F'\t' '{n=split($1,a,"/"); name=a[n]; sub(/\.fna$/,"",name); print name"\t"$6}' \
    | sort | uniq -c | sort -rn | head -20 | tee -a "$LOG" || true
