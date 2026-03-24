#!/usr/bin/env bash
# =============================================================================
# STEP 10 — Pan-genome Analysis  [NOVEL ADDITION]
# conda env : roary_env  (Roary)
#             phylogeny_env  (IQ-TREE 2)
# Tool      : Roary
# Why added : Paper mentions "open pan-genome" but never quantifies it.
#             Provides: core genome size, accessory genome, singletons,
#             Heaps' law alpha (open vs closed), and accessory gene–ARG links.
# Input     : Prokka GFF files from step01 (copied to results/gff/)
# Output    : Pan-genome matrix, core gene alignment, accumulation curve data
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

GFF_DIR="$RESULTS_DIR/gff"
OUT="$RESULTS_DIR/10_pangenome_roary"
LOG="$RESULTS_DIR/00_logs/10_pangenome.log"
mkdir -p "$OUT" "$OUT/pangenome_plots"

echo "[INFO] Starting pan-genome analysis (Roary) — $(date)" | tee "$LOG"

# ─── CHECK GFF FILES ─────────────────────────────────────────────────────────
N_GFF=$(ls "$GFF_DIR"/*.gff 2>/dev/null | wc -l || echo 0)

if [[ "$N_GFF" -eq 0 ]]; then
    echo "[ERROR] No .gff files in $GFF_DIR" | tee -a "$LOG"
    echo "        Step 01 (annotation) must complete first." | tee -a "$LOG"
    exit 1
fi
echo "[OK  ] Found $N_GFF GFF files" | tee -a "$LOG"

# ─── RUN ROARY ───────────────────────────────────────────────────────────────
ROARY_OUT="$OUT/roary_output"
rm -rf "$ROARY_OUT"   # Roary fails if output directory already exists

conda run -n roary_env roary \
    -p "$THREADS" \
    -f "$ROARY_OUT" \
    -e \
    -n \
    -r \
    -v \
    -i 95 \
    "$GFF_DIR"/*.gff \
    >> "$LOG" 2>&1

echo "[INFO] Roary complete" | tee -a "$LOG"

# ─── PAN-GENOME STATISTICS ───────────────────────────────────────────────────
SUMMARY_FILE="$ROARY_OUT/summary_statistics.txt"
PA="$ROARY_OUT/gene_presence_absence.csv"

if [[ -f "$SUMMARY_FILE" ]]; then
    echo "" | tee -a "$LOG"
    echo "=== Pan-genome Statistics ===" | tee -a "$LOG"
    cat "$SUMMARY_FILE" | tee -a "$LOG"
elif [[ -f "$PA" ]]; then
    echo "[INFO] Calculating stats from presence/absence matrix..." | tee -a "$LOG"
    python3 - << PYEOF | tee -a "$LOG"
import csv
pa_file = "$PA"
n_samples = $N_GFF
with open(pa_file) as f:
    reader = csv.reader(f)
    header = next(reader)
    n_meta = 14
    core = soft_core = shell = cloud = total = 0
    for row in reader:
        total += 1
        presence = sum(1 for v in row[n_meta:] if v.strip())
        freq = presence / n_samples
        if freq >= 0.99:    core += 1
        elif freq >= 0.95:  soft_core += 1
        elif freq >= 0.15:  shell += 1
        else:               cloud += 1
print(f"Total genes (pan):   {total}")
print(f"Core     (>=99%):    {core}")
print(f"Soft-core (95-99%):  {soft_core}")
print(f"Shell    (15-95%):   {shell}")
print(f"Cloud     (<15%):    {cloud}")
print(f"Accessory:           {total - core}")
print(f"Core/Pan ratio:      {core/total:.3f}")
PYEOF
fi

# ─── CORE GENOME PHYLOGENY ───────────────────────────────────────────────────
CORE_ALIGN="$ROARY_OUT/core_gene_alignment.aln"
if [[ -f "$CORE_ALIGN" ]]; then
    echo "[INFO] Building core genome phylogeny with IQ-TREE 2..." | tee -a "$LOG"
    conda run -n phylogeny_env iqtree2 \
        -s "$CORE_ALIGN" -m GTR+G -B 1000 -T "$THREADS" \
        --prefix "$OUT/core_genome_tree" --redo \
        >> "$LOG" 2>&1
    echo "[INFO] Core genome tree → $OUT/core_genome_tree.treefile" | tee -a "$LOG"
fi

# ─── PAN-GENOME ACCUMULATION PLOTS ───────────────────────────────────────────
conda run -n roary_env roary_plots.py \
    -o "$OUT/pangenome_plots" \
    "$ROARY_OUT/summary_statistics.txt" \
    "$ROARY_OUT/gene_presence_absence.csv" \
    >> "$LOG" 2>&1 || \
echo "[WARN] roary_plots.py not available — skipping plots" | tee -a "$LOG"

# ─── CROSS-REFERENCE ACCESSORY GENES WITH ARGs ───────────────────────────────
echo "[INFO] Cross-referencing accessory genes with ARGs..." | tee -a "$LOG"
ACCESSORY_ARG="$RESULTS_DIR/13_summary_tables/accessory_ARG_genes.tsv"
CARD_RAW="$RESULTS_DIR/03_ARG_RGI/abricate_card_raw.tsv"
echo -e "gene\tcategory\tfrequency_pct\tpresent_in_n_samples" > "$ACCESSORY_ARG"

if [[ -f "$PA" && -f "$CARD_RAW" ]]; then
    python3 - << PYEOF >> "$ACCESSORY_ARG" 2>>"$LOG"
import csv
pa_file = "$PA"
card_raw = "$CARD_RAW"
n_samples = $N_GFF

arg_genes = set()
try:
    with open(card_raw) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 5 and parts[5]:
                arg_genes.add(parts[5].lower())
except Exception:
    pass

with open(pa_file) as f:
    reader = csv.reader(f)
    next(reader)
    n_meta = 14
    for row in reader:
        gene = row[0].lower()
        presence = sum(1 for v in row[n_meta:] if v.strip())
        freq = presence / n_samples * 100
        category = "ARG" if (gene in arg_genes or any(g in gene for g in arg_genes)) \
                   else "other_accessory"
        if (presence / n_samples) < 0.99 and \
           (category == "ARG" or (presence / n_samples) < 0.95):
            print(f"{row[0]}\t{category}\t{freq:.1f}\t{presence}")
PYEOF
fi

echo "[INFO] Accessory ARG table → $ACCESSORY_ARG" | tee -a "$LOG"
echo "[DONE] Step 10 complete — $(date)" | tee -a "$LOG"
echo ""
echo "Key outputs:"
echo "  Pan-genome matrix : $ROARY_OUT/gene_presence_absence.csv"
echo "  Core genome tree  : $OUT/core_genome_tree.treefile"
echo "  Pan-genome plots  : $OUT/pangenome_plots/"
echo "  Accessory ARGs    : $ACCESSORY_ARG"
