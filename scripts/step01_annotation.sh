#!/usr/bin/env bash
# =============================================================================
# STEP 01 — Genome Annotation with Prokka
# conda env : annotation_env
# Tool      : prokka 1.14.x
# Paper ref : "Annotations carried out using Prokka v2.1.1"
# Output    : GFF, GBK, FAA, FFN per genome
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

OUT="$RESULTS_DIR/01_annotation"
LOG="$RESULTS_DIR/00_logs/01_annotation.log"
mkdir -p "$OUT"

GENUS="Acinetobacter"
SPECIES="baumannii"

echo "[INFO] Starting Prokka annotation — $(date)" | tee "$LOG"
echo "[INFO] Threads: $THREADS" | tee -a "$LOG"

while IFS= read -r SAMPLE; do
    GENOME="$ASSEMBLY_DIR/${SAMPLE}.fna"
    OUTDIR="$OUT/$SAMPLE"

    if [[ ! -f "$GENOME" ]]; then
        echo "[WARN] Genome not found: $GENOME — skipping" | tee -a "$LOG"
        continue
    fi

    if [[ -f "$OUTDIR/${SAMPLE}.gff" ]]; then
        echo "[SKIP] $SAMPLE already annotated" | tee -a "$LOG"
        continue
    fi

    echo "[RUN ] Annotating: $SAMPLE" | tee -a "$LOG"

    conda run -n annotation_env prokka \
        --outdir "$OUTDIR" \
        --prefix "$SAMPLE" \
        --genus "$GENUS" \
        --species "$SPECIES" \
        --strain "$SAMPLE" \
        --kingdom Bacteria \
        --cpus "$THREADS" \
        --compliant \
        --force \
        "$GENOME" >> "$LOG" 2>&1

    echo "[DONE] $SAMPLE" | tee -a "$LOG"
done < "$SAMPLE_LIST"

# ─── COLLECT SUMMARY TABLE ────────────────────────────────────────────────────
echo "[INFO] Generating annotation summary table..."
SUMMARY="$RESULTS_DIR/13_summary_tables/annotation_summary.tsv"
echo -e "sample\tCDS\ttRNA\trRNA\tgenome_size_bp\tGC_pct" > "$SUMMARY"

while IFS= read -r SAMPLE; do
    TXT="$OUT/$SAMPLE/${SAMPLE}.txt"
    [[ ! -f "$TXT" ]] && continue

    CDS=$(grep  "^CDS:"  "$TXT" | awk '{print $2}')
    TRNA=$(grep "^tRNA:" "$TXT" | awk '{print $2}')
    RRNA=$(grep "^rRNA:" "$TXT" | awk '{print $2}')
    GENOME="$ASSEMBLY_DIR/${SAMPLE}.fna"
    SIZE=$(grep -v "^>" "$GENOME" | tr -d '\n' | wc -c)
    GC=$(grep -v "^>" "$GENOME" | tr -d '\n' | \
         awk '{c=gsub(/[GCgc]/,"",$0); print c/length($0)*100}')

    echo -e "${SAMPLE}\t${CDS}\t${TRNA}\t${RRNA}\t${SIZE}\t${GC}" >> "$SUMMARY"
done < "$SAMPLE_LIST"

# ─── COPY GFF FILES TO SHARED LOCATION FOR ROARY (Step 10) ──────────────────
GFF_DIR="$RESULTS_DIR/gff"
mkdir -p "$GFF_DIR"
while IFS= read -r SAMPLE; do
    GFF="$OUT/$SAMPLE/${SAMPLE}.gff"
    [[ -f "$GFF" ]] && cp "$GFF" "$GFF_DIR/${SAMPLE}.gff"
done < "$SAMPLE_LIST"
echo "[INFO] GFF files copied to: $GFF_DIR"

echo "[INFO] Annotation summary → $SUMMARY"
echo "[DONE] Step 01 complete — $(date)" | tee -a "$LOG"
