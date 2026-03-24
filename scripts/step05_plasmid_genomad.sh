#!/usr/bin/env bash
# =============================================================================
# STEP 05 — Plasmid Detection and Typing
# conda env : clinical_env  (geNomad, MOB-suite, abricate)
# Tools     : geNomad 1.8.x  — neural-network plasmid/virus classification
#             MOB-suite       — replicon typing
# Paper ref : PlasmidFinder / replicon typing
# Output    : Per-sample plasmid FASTA, summary TSV, MOB-suite replicon table
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

OUT="$RESULTS_DIR/06_plasmid_genomad"
LOG="$RESULTS_DIR/00_logs/05_plasmid.log"
mkdir -p "$OUT"

echo "[INFO] Starting plasmid detection — $(date)" | tee "$LOG"
echo "[INFO] geNomad DB: $GENOMAD_DB" | tee -a "$LOG"

if [[ ! -d "$GENOMAD_DB" ]]; then
    echo "[ERROR] geNomad database not found: $GENOMAD_DB" | tee -a "$LOG"
    echo "        Run: genomad download-database <dir>" | tee -a "$LOG"
    exit 1
fi

# ─── VERIFY mmseqs2 ──────────────────────────────────────────────────────────
if ! conda run -n clinical_env mmseqs version &>/dev/null 2>&1; then
    echo "[INST] Installing mmseqs2 into clinical_env..." | tee -a "$LOG"
    conda install -n clinical_env -c conda-forge -c bioconda mmseqs2 -y >> "$LOG" 2>&1
fi

# ─── CLEANUP INCOMPLETE RUNS ─────────────────────────────────────────────────
while IFS= read -r SAMPLE; do
    SAMPLE_OUT="$OUT/$SAMPLE"
    DONE_MARKER="$SAMPLE_OUT/${SAMPLE}_summary/${SAMPLE}_plasmid_summary.tsv"
    if [[ -d "$SAMPLE_OUT/${SAMPLE}_annotate" && ! -f "$DONE_MARKER" ]]; then
        rm -rf "$SAMPLE_OUT"
    fi
done < "$SAMPLE_LIST"

# ─── MAIN LOOP ───────────────────────────────────────────────────────────────
PLASMID_SUMMARY="$RESULTS_DIR/13_summary_tables/plasmid_summary.tsv"
echo -e "sample\tn_plasmid_contigs\tn_plasmid_ARGs\ttotal_plasmid_bp" > "$PLASMID_SUMMARY"

while IFS= read -r SAMPLE; do
    GENOME="$ASSEMBLY_DIR/${SAMPLE}.fna"
    SAMPLE_OUT="$OUT/$SAMPLE"
    DONE_MARKER="$SAMPLE_OUT/${SAMPLE}_summary/${SAMPLE}_plasmid_summary.tsv"

    [[ ! -f "$GENOME" ]] && continue

    if [[ -f "$DONE_MARKER" ]]; then
        echo "[SKIP] $SAMPLE already processed." | tee -a "$LOG"
    else
        echo "[RUN ] geNomad: $SAMPLE" | tee -a "$LOG"
        mkdir -p "$SAMPLE_OUT"

        if ! conda run -n clinical_env genomad end-to-end \
                --cleanup --splits 8 --threads "$THREADS" \
                --min-score 0.7 --enable-score-calibration \
                "$GENOME" "$SAMPLE_OUT" "$GENOMAD_DB" >> "$LOG" 2>&1; then
            echo "[WARN] geNomad failed for $SAMPLE." | tee -a "$LOG"
            continue
        fi
    fi

    PLASMID_FA="$SAMPLE_OUT/${SAMPLE}_summary/${SAMPLE}_plasmid.fna"
    N_CONTIGS=0; TOTAL_BP=0; N_ARG=0

    if [[ -f "$PLASMID_FA" && -s "$PLASMID_FA" ]]; then
        N_CONTIGS=$(grep -c "^>" "$PLASMID_FA" || echo 0)
        TOTAL_BP=$(grep -v "^>" "$PLASMID_FA" | tr -d '\n' | wc -c || echo 0)
        cp "$PLASMID_FA" "$OUT/${SAMPLE}_plasmids.fna"

        if [[ "$N_CONTIGS" -gt 0 ]]; then
            conda run -n clinical_env abricate \
                --db card --minid 80 --mincov 60 --threads "$THREADS" \
                "$OUT/${SAMPLE}_plasmids.fna" \
                > "$OUT/${SAMPLE}_plasmid_ARGs.tsv" 2>>"$LOG" || true
            N_ARG=$(awk 'NR>1' "$OUT/${SAMPLE}_plasmid_ARGs.tsv" 2>/dev/null | wc -l || echo 0)
        fi
    fi

    echo -e "${SAMPLE}\t${N_CONTIGS}\t${N_ARG}\t${TOTAL_BP}" >> "$PLASMID_SUMMARY"
    echo "[DONE] $SAMPLE — ${N_CONTIGS} plasmids, ${N_ARG} ARGs" | tee -a "$LOG"
done < "$SAMPLE_LIST"

# ─── MOB-SUITE REPLICON TYPING ───────────────────────────────────────────────
echo "[INFO] Running MOB-suite replicon typing..." | tee -a "$LOG"
MOBOUT="$RESULTS_DIR/11_mob_suite"
mkdir -p "$MOBOUT"

while IFS= read -r SAMPLE; do
    GENOME="$ASSEMBLY_DIR/${SAMPLE}.fna"
    [[ ! -f "$GENOME" ]] && continue
    [[ -f "$MOBOUT/${SAMPLE}_mobtyper.tsv" ]] && continue

    echo "[RUN ] MOB-suite: $SAMPLE" | tee -a "$LOG"
    conda run -n clinical_env mob_typer \
        --infile "$GENOME" \
        --out_file "$MOBOUT/${SAMPLE}_mobtyper.tsv" \
        --num_threads "$THREADS" \
        >> "$LOG" 2>&1 || \
    echo "[WARN] MOB-suite failed: $SAMPLE" | tee -a "$LOG"
done < "$SAMPLE_LIST"

# ─── AGGREGATE MOB-SUITE RESULTS ─────────────────────────────────────────────
MOB_MERGED="$RESULTS_DIR/13_summary_tables/mob_suite_summary.tsv"
FIRST_MOB=$(ls "$MOBOUT"/*_mobtyper.tsv 2>/dev/null | head -1 || true)
if [[ -n "$FIRST_MOB" ]]; then
    head -1 "$FIRST_MOB" > "$MOB_MERGED"
    for f in "$MOBOUT"/*_mobtyper.tsv; do
        [[ -f "$f" ]] && tail -n +2 "$f" >> "$MOB_MERGED"
    done
fi

cat "$OUT"/*_plasmids.fna > "$OUT/all_plasmids_merged.fna" 2>/dev/null || true

echo "" | tee -a "$LOG"
echo "=== Plasmid Summary ===" | tee -a "$LOG"
column -t "$PLASMID_SUMMARY" | tee -a "$LOG"
echo "[DONE] Step 05 complete — $(date)" | tee -a "$LOG"
