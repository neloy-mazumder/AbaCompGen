#!/usr/bin/env bash
# =============================================================================
# STEP 11 — GTDB-Tk Phylogenomics + Temporal ARG Trend Analysis  [NOVEL]
# conda env : assembly_binning_sn_env  (gtdbtk 2.6.1)
# Tool      : GTDB-Tk classify_wf
# Why added : Provides GTDB-based taxonomic placement and reference-anchored
#             phylogenomic tree — more robust than SNP-only for deep diversity.
#             Temporal ARG analysis (acquisition over time) extends the paper.
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

GTDB_OUT="$RESULTS_DIR/12_gtdbtk"
TEMPORAL_OUT="$RESULTS_DIR/13_summary_tables"
LOG="$RESULTS_DIR/00_logs/11_gtdbtk_temporal.log"
mkdir -p "$GTDB_OUT"

echo "[INFO] Starting GTDB-Tk + temporal analysis — $(date)" | tee "$LOG"

# ─── PART A: GTDB-Tk Classify ────────────────────────────────────────────────
echo "[INFO] === PART A: GTDB-Tk classify ===" | tee -a "$LOG"
echo "[INFO] GTDB data path: $GTDB_DATA" | tee -a "$LOG"

if [[ ! -d "$GTDB_DATA" ]]; then
    echo "[ERROR] GTDB database not found: $GTDB_DATA" | tee -a "$LOG"
    echo "        Download from https://gtdb.ecogenomic.org/downloads" | tee -a "$LOG"
    exit 1
fi

export GTDBTK_DATA_PATH="$GTDB_DATA"

conda run -n assembly_binning_sn_env gtdbtk classify_wf \
    --genome_dir "$ASSEMBLY_DIR" \
    --extension fna \
    --out_dir "$GTDB_OUT" \
    --cpus "$THREADS" \
    --skip_ani_screen \
    >> "$LOG" 2>&1 || \
echo "[WARN] GTDB-Tk failed — check database path: $GTDB_DATA" | tee -a "$LOG"

if [[ -f "$GTDB_OUT/gtdbtk.bac120.summary.tsv" ]]; then
    echo "" | tee -a "$LOG"
    echo "=== GTDB taxonomy (first 5 samples) ===" | tee -a "$LOG"
    awk -F'\t' 'NR<=6{print $1"\t"$2}' \
        "$GTDB_OUT/gtdbtk.bac120.summary.tsv" | tee -a "$LOG"
fi

# ─── PART B: Temporal ARG Accumulation Analysis ───────────────────────────────
echo "" | tee -a "$LOG"
echo "[INFO] === PART B: Temporal ARG trend analysis ===" | tee -a "$LOG"

TEMPORAL_TABLE="$TEMPORAL_OUT/temporal_ARG_trends.tsv"
echo -e "sample\tyear\thospital\tn_total_ARG\tn_acquired_ARG\tn_intrinsic_ARG\tn_bla_genes\tn_aminoglycoside\tn_sul_genes\thas_blaOXA23" \
    > "$TEMPORAL_TABLE"

# A. baumannii intrinsic ARGs (always present regardless of resistance phenotype)
INTRINSIC_ARGS="adeABC|adeIJK|adeR|adeS|oxa-51|blaOXA-51|ADC|blaADC|gyrA|parC|lpsB"
CARD_RAW="$RESULTS_DIR/03_ARG_RGI/abricate_card_raw.tsv"

while IFS= read -r SAMPLE; do
    YEAR=$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $2}' "$METADATA" 2>/dev/null || echo "NA")
    HOSPITAL=$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $3}' "$METADATA" 2>/dev/null || echo "NA")

    N_TOTAL=0; N_ACQUIRED=0; N_INTRINSIC=0
    N_BLA=0; N_AMG=0; N_SUL=0; HAS_OXA23="no"

    if [[ -f "$CARD_RAW" ]]; then
        SAMPLE_LINES=$(grep "/${SAMPLE}\.fna" "$CARD_RAW" 2>/dev/null || true)
        N_TOTAL=$(echo "$SAMPLE_LINES" | grep -c "." || echo 0)
        N_INTRINSIC=$(echo "$SAMPLE_LINES" | grep -ciE "$INTRINSIC_ARGS" || echo 0)
        N_ACQUIRED=$(( N_TOTAL > N_INTRINSIC ? N_TOTAL - N_INTRINSIC : 0 ))
        N_BLA=$(echo "$SAMPLE_LINES" | grep -ciE "bla|OXA|ADC|TEM|NDM|IMP|VIM" || echo 0)
        N_AMG=$(echo "$SAMPLE_LINES" | grep -ciE "aac|aph|ant|aad|arm" || echo 0)
        N_SUL=$(echo "$SAMPLE_LINES" | grep -ciE "sul1|sul2" || echo 0)
        echo "$SAMPLE_LINES" | grep -iq "OXA-23\|OXA23\|blaOXA-23" && HAS_OXA23="yes" || true
    fi

    echo -e "${SAMPLE}\t${YEAR}\t${HOSPITAL}\t${N_TOTAL}\t${N_ACQUIRED}\t${N_INTRINSIC}\t${N_BLA}\t${N_AMG}\t${N_SUL}\t${HAS_OXA23}" \
        >> "$TEMPORAL_TABLE"
done < "$SAMPLE_LIST"

echo "[INFO] Temporal ARG table → $TEMPORAL_TABLE" | tee -a "$LOG"

# ─── TEMPORAL TREND SUMMARY ──────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "=== ARG Trends by Year ===" | tee -a "$LOG"
awk -F'\t' 'NR>1 && $2 != "NA" {
    year[$2] += $5; count[$2]++
} END {
    print "Year\tMean_Acquired_ARGs\tN_genomes"
    for (y in year) printf "%s\t%.2f\t%d\n", y, year[y]/count[y], count[y]
}' "$TEMPORAL_TABLE" | sort -k1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== blaOXA-23 Prevalence by Year ===" | tee -a "$LOG"
awk -F'\t' 'NR>1 && $2 != "NA" {
    total[$2]++
    if ($10=="yes") oxa23[$2]++
} END {
    print "Year\tTotal\tWith_blaOXA23\tPct"
    for (y in total) printf "%s\t%d\t%d\t%.1f%%\n", y, total[y], oxa23[y]+0, (oxa23[y]+0)/total[y]*100
}' "$TEMPORAL_TABLE" | sort -k1 | tee -a "$LOG"

# ─── PART C: Hospital-level ARG burden ───────────────────────────────────────
echo "" | tee -a "$LOG"
echo "=== Hospital-level ARG Burden ===" | tee -a "$LOG"
awk -F'\t' 'NR>1 && $3 != "NA" {
    hosp[$3] += $5; count[$3]++
} END {
    print "Hospital\tMean_Acquired_ARGs\tN_genomes"
    for (h in hosp) printf "%s\t%.2f\t%d\n", h, hosp[h]/count[h], count[h]
}' "$TEMPORAL_TABLE" | sort -k2 -rn | tee -a "$LOG"

echo "[DONE] Step 11 complete — $(date)" | tee -a "$LOG"
