#!/usr/bin/env bash
# =============================================================================
# STEP 06 — Insertion Sequence (IS) Detection
# conda env : clinical_env
# Tool      : ISEScan
# Paper ref : "BLAST against ISfinder — 16 IS families detected"
# Bonus     : IS-ARG linkage analysis (IS elements within 5 kb of ARGs)
# Output    : Per-sample IS table + cross-sample IS family matrix
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

OUT="$RESULTS_DIR/07_IS_isescan"
LOG="$RESULTS_DIR/00_logs/06_IS.log"
mkdir -p "$OUT"

echo "[INFO] Starting IS detection (ISEScan) — $(date)" | tee "$LOG"

# ─── PATCH: fastcluster.py for numpy >= 1.24 ─────────────────────────────────
# ISEScan calls fastcluster.linkage(..., preserve_input=False) which triggers
# numpy.array(X, copy=None) — rejected in numpy >= 1.24. Fix: copy=False.
FASTCLUSTER_PY=$(conda run -n clinical_env python -c \
    "import fastcluster, os; print(os.path.abspath(fastcluster.__file__))" \
    2>/dev/null || echo "")

if [[ -n "$FASTCLUSTER_PY" && -f "$FASTCLUSTER_PY" ]]; then
    if grep -q "copy=True if preserve_input else None" "$FASTCLUSTER_PY"; then
        echo "[FIX ] Patching fastcluster.py (numpy >= 1.24 compatibility)..." | tee -a "$LOG"
        sed -i 's/copy=True if preserve_input else None/copy=True if preserve_input else False/g' \
            "$FASTCLUSTER_PY"
        echo "[OK  ] fastcluster.py patched." | tee -a "$LOG"
    else
        echo "[OK  ] fastcluster.py already compatible." | tee -a "$LOG"
    fi
else
    echo "[WARN] Could not locate fastcluster.py — ISEScan may fail." | tee -a "$LOG"
fi

# ─── RUN ISEScan ─────────────────────────────────────────────────────────────
while IFS= read -r SAMPLE; do
    GENOME="$ASSEMBLY_DIR/${SAMPLE}.fna"
    SAMPLE_OUT="$OUT/$SAMPLE"
    [[ ! -f "$GENOME" ]] && continue

    if [[ -f "$SAMPLE_OUT/${SAMPLE}.fna.tsv" ]]; then
        echo "[SKIP] $SAMPLE already processed" | tee -a "$LOG"
        continue
    fi

    echo "[RUN ] ISEScan: $SAMPLE" | tee -a "$LOG"
    mkdir -p "$SAMPLE_OUT"
    cp "$GENOME" "$SAMPLE_OUT/${SAMPLE}.fna"

    conda run -n clinical_env isescan.py \
        --seqfile "$SAMPLE_OUT/${SAMPLE}.fna" \
        --output "$SAMPLE_OUT" \
        --nthread "$THREADS" \
        >> "$LOG" 2>&1 || \
    echo "[WARN] ISEScan failed for $SAMPLE" | tee -a "$LOG"

    echo "[DONE] $SAMPLE" | tee -a "$LOG"
done < "$SAMPLE_LIST"

# ─── MERGE ALL IS RESULTS ─────────────────────────────────────────────────────
echo "[INFO] Merging IS results..." | tee -a "$LOG"
IS_MERGED="$OUT/all_IS_elements.tsv"
FIRST_TSV=$(find "$OUT" -name "*.tsv" | grep -v "all_IS" | head -1 || true)

if [[ -f "$FIRST_TSV" ]]; then
    head -1 "$FIRST_TSV" | sed 's/^/sample\t/' > "$IS_MERGED"
    while IFS= read -r SAMPLE; do
        TSV="$OUT/$SAMPLE/${SAMPLE}.fna.tsv"
        [[ ! -f "$TSV" ]] && continue
        tail -n +2 "$TSV" | awk -v s="$SAMPLE" '{print s"\t"$0}' >> "$IS_MERGED"
    done < "$SAMPLE_LIST"
fi

# ─── IS FAMILY FREQUENCY TABLE ───────────────────────────────────────────────
IS_FREQ="$RESULTS_DIR/13_summary_tables/IS_family_frequency.tsv"
echo -e "sample\tIS_family\tcount" > "$IS_FREQ"

while IFS= read -r SAMPLE; do
    SAMPLE=$(echo "$SAMPLE" | tr -d '\r' | xargs)
    TSV="$OUT/$SAMPLE/${SAMPLE}.fna.tsv"
    [[ ! -f "$TSV" ]] && TSV="$OUT/$SAMPLE/$SAMPLE/${SAMPLE}.fna.tsv"
    [[ ! -f "$TSV" ]] && continue

    tail -n +2 "$TSV" | awk -v s="$SAMPLE" '
    $2 != "" { fam[$2]++ }
    END { for (f in fam) print s"\t"f"\t"fam[f] }' >> "$IS_FREQ"
done < "$SAMPLE_LIST"

echo "[INFO] IS family frequency → $IS_FREQ" | tee -a "$LOG"

# ─── IS-ARG LINKAGE (IS within 5 kb of an ARG) ───────────────────────────────
IS_ARG_LINK="$RESULTS_DIR/13_summary_tables/IS_ARG_linkage.tsv"
echo -e "sample\tIS_family\tIS_start\tIS_end\tARG\tARG_start\tARG_end\tdistance_bp" > "$IS_ARG_LINK"

CARD_RAW="$RESULTS_DIR/03_ARG_RGI/abricate_card_raw.tsv"

while IFS= read -r SAMPLE; do
    SAMPLE=$(echo "$SAMPLE" | tr -d '\r' | xargs)
    IS_TSV="$OUT/$SAMPLE/${SAMPLE}.fna.tsv"
    [[ ! -f "$IS_TSV" ]] && IS_TSV="$OUT/$SAMPLE/$SAMPLE/${SAMPLE}.fna.tsv"
    [[ ! -f "$IS_TSV" ]] && continue

    awk -v cur_sample="$SAMPLE" -v card_raw="$CARD_RAW" '
    BEGIN {
        idx = 0
        while ((getline line < card_raw) > 0) {
            if (line ~ /^#/) continue
            split(line, a, "\t")
            n = split(a[1], p, "/"); fname = p[n]
            sub(/\.fna$/, "", fname); sub(/\.fasta$/, "", fname)
            if (fname == cur_sample) {
                ac = a[2]; gsub(/[[:space:]]/, "", ac)
                arg_contig[idx] = ac
                arg_start[idx]  = a[3]
                arg_end[idx]    = a[4]
                arg_gene[idx]   = a[6]
                idx++
            }
        }
        close(card_raw)
    }
    {
        ic = $1; gsub(/[[:space:]]/, "", ic)
        is_fam = $2; is_start = $4; is_end = $5
        for (i=0; i<idx; i++) {
            if (ic == arg_contig[i]) {
                d1 = is_start - arg_end[i]
                d2 = arg_start[i] - is_end
                dist = (d1 > 0) ? d1 : (d2 > 0 ? d2 : 0)
                if (dist <= 5000) {
                    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n", \
                    cur_sample, is_fam, is_start, is_end, \
                    arg_gene[i], arg_start[i], arg_end[i], dist
                }
            }
        }
    }' "$IS_TSV" >> "$IS_ARG_LINK"
done < "$SAMPLE_LIST"

echo "[INFO] IS-ARG linkage table → $IS_ARG_LINK" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== IS Family Distribution ===" | tee -a "$LOG"
awk -F'\t' 'NR>1{print $2}' "$IS_FREQ" \
    | sort | uniq -c | sort -rn | head -20 | tee -a "$LOG" || true

echo "" | tee -a "$LOG"
echo "=== Most Frequent IS-ARG Linkages ===" | tee -a "$LOG"
awk -F'\t' 'NR>1{print $2"\t"$5}' "$IS_ARG_LINK" \
    | sort | uniq -c | sort -rn | head -20 | tee -a "$LOG" || true

echo "[DONE] Step 06 complete — $(date)" | tee -a "$LOG"
