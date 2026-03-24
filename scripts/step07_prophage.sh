#!/usr/bin/env bash
# =============================================================================
# STEP 07 — Prophage Detection
# conda env : clinical_env  (geNomad viral mode)
# Tool      : geNomad end-to-end (viral classification)
# Paper ref : "PHASTER — 274 prophage regions"
# Bonus     : Checks whether prophage regions carry ARGs
#             Optional PHASTER API submission script included
# Output    : Per-sample prophage FASTA, summary TSV, ARG-in-phage check
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

OUT="$RESULTS_DIR/08_prophage_genomad"
LOG="$RESULTS_DIR/00_logs/07_prophage.log"
PHAGE_SUMMARY="$RESULTS_DIR/13_summary_tables/prophage_summary.tsv"
mkdir -p "$OUT" "$RESULTS_DIR/13_summary_tables"

echo "[INFO] Running geNomad prophage detection — $(date)" | tee "$LOG"
echo "[INFO] geNomad DB: $GENOMAD_DB" | tee -a "$LOG"

if [[ ! -d "$GENOMAD_DB" ]]; then
    echo "[ERROR] geNomad DB not found: $GENOMAD_DB" | tee -a "$LOG"
    echo "        Run: genomad download-database <dir>" | tee -a "$LOG"
    exit 1
fi

echo -e "sample\tn_virus_contigs\tn_intact\tn_provirus\ttotal_phage_bp" > "$PHAGE_SUMMARY"

while IFS= read -r SAMPLE; do
    ASM="$ASSEMBLY_DIR/${SAMPLE}.fna"
    if [[ ! -f "$ASM" ]]; then
        echo "[WARN] $SAMPLE assembly missing" | tee -a "$LOG"
        echo -e "${SAMPLE}\t0\t0\t0\t0" >> "$PHAGE_SUMMARY"
        continue
    fi

    GENOMAD_OUTDIR="$OUT/${SAMPLE}_genomad"
    VSUM="$GENOMAD_OUTDIR/${SAMPLE}_summary/${SAMPLE}_virus_summary.tsv"
    VFASTA="$GENOMAD_OUTDIR/${SAMPLE}_summary/${SAMPLE}_virus.fna"

    if [[ ! -f "$VSUM" ]]; then
        echo "[RUN ] geNomad: $SAMPLE" | tee -a "$LOG"
        conda run -n clinical_env genomad end-to-end \
            --cleanup --splits 8 \
            "$ASM" "$GENOMAD_OUTDIR" "$GENOMAD_DB" \
            >> "$LOG" 2>&1
    else
        echo "[SKIP] $SAMPLE already processed." | tee -a "$LOG"
    fi

    if [[ -f "$VSUM" ]]; then
        N_VIRUS=$(awk 'NR>1 {print}' "$VSUM" | wc -l)
        N_PRO=$(awk 'NR>1 && /provirus/ {print}' "$VSUM" | wc -l)
        TOTAL_BP=$(awk 'NR>1 {sum+=$2} END {print sum+0}' "$VSUM")
        cp "$VFASTA" "$OUT/${SAMPLE}_prophage.fna" 2>/dev/null || true
    else
        N_VIRUS=0; N_PRO=0; TOTAL_BP=0
        touch "$OUT/${SAMPLE}_prophage.fna"
    fi

    echo -e "${SAMPLE}\t${N_VIRUS}\t0\t${N_PRO}\t${TOTAL_BP}" >> "$PHAGE_SUMMARY"
    echo "[INFO] $SAMPLE — $N_VIRUS virus contigs ($TOTAL_BP bp)" | tee -a "$LOG"
done < "$SAMPLE_LIST"

# ─── ARGs IN PROPHAGE REGIONS ────────────────────────────────────────────────
echo "[INFO] Checking for ARGs within prophage regions..." | tee -a "$LOG"
PHAGE_ARG="$RESULTS_DIR/13_summary_tables/prophage_ARG_check.tsv"
echo -e "sample\tphage_contig\tARG_gene\tpident\tcoverage" > "$PHAGE_ARG"

while IFS= read -r SAMPLE; do
    PHAGE_FA="$OUT/${SAMPLE}_prophage.fna"
    [[ ! -f "$PHAGE_FA" || $(wc -c < "$PHAGE_FA") -lt 100 ]] && continue
    conda run -n clinical_env abricate \
        --db card --minid 80 --mincov 60 --threads "$THREADS" \
        "$PHAGE_FA" 2>>"$LOG" | \
    awk -F'\t' -v s="$SAMPLE" 'NR>1{print s"\t"$2"\t"$6"\t"$10"\t"$9}' \
        >> "$PHAGE_ARG" || true
done < "$SAMPLE_LIST"

N_PHAGE_ARG=$(awk 'NR>1' "$PHAGE_ARG" | wc -l)
echo "[INFO] ARGs in prophage regions: $N_PHAGE_ARG" | tee -a "$LOG"

# ─── OPTIONAL: PHASTER API SUBMISSION ────────────────────────────────────────
PHASTER_SCRIPT="$OUT/submit_to_phaster_api.sh"
cat > "$PHASTER_SCRIPT" << 'PHEOF'
#!/usr/bin/env bash
# Optional: Submit assemblies to PHASTER web API for intact prophage analysis.
# Requires internet access. Cross-validates with the paper's exact tool.
# API docs: https://phaster.ca/
set -euo pipefail
source "$PROJECT_DIR/env_config.sh"
PHASTER_OUT="$RESULTS_DIR/08_prophage_genomad/phaster_api"
mkdir -p "$PHASTER_OUT"
for GENOME in "$ASSEMBLY_DIR"/*.fna; do
    SAMPLE=$(basename "$GENOME" .fna)
    echo "[SUBMIT] $SAMPLE to PHASTER API..."
    JOB_ID=$(curl -s -X POST --form "contigs=@${GENOME}" \
        "https://phaster.ca/phaster_api?acc=${SAMPLE}" \
        | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('job_id',''))")
    echo "$SAMPLE $JOB_ID" >> "$PHASTER_OUT/job_ids.txt"
    echo "[INFO] $SAMPLE → Job ID: $JOB_ID"
    sleep 2
done
echo "Retrieve results: https://phaster.ca/"
PHEOF
chmod +x "$PHASTER_SCRIPT"
echo "[INFO] PHASTER submission script → $PHASTER_SCRIPT (optional, requires internet)"

echo "" | tee -a "$LOG"
echo "=== Prophage Summary ===" | tee -a "$LOG"
awk -F'\t' 'NR>1{total+=$5; nv+=$2; np+=$4}
END{print "Total virus contigs: "nv"\nProviruses: "np"\nTotal phage bp: "total}' \
    "$PHAGE_SUMMARY" | tee -a "$LOG"

echo "[DONE] Step 07 complete — $(date)" | tee -a "$LOG"
