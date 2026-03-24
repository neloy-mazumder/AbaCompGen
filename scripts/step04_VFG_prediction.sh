#!/usr/bin/env bash
# =============================================================================
# STEP 04 — Virulence Factor Gene (VFG) Prediction
# conda env : clinical_env
# Tools     : abricate (VFDB)
#             BLAST against VFDB core dataset (mirrors VFanalyzer)
# Paper ref : "VFanalyzer against VFDB — 56 VFG types in 9 categories"
# Output    : Per-sample + merged presence/absence matrix + category summary
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

OUT="$RESULTS_DIR/05_VFG_abricate"
LOG="$RESULTS_DIR/00_logs/04_VFG.log"
mkdir -p "$OUT"

echo "[INFO] Starting VFG prediction — $(date)" | tee "$LOG"

# ─── PART A: Abricate + VFDB ─────────────────────────────────────────────────
echo "[INFO] === PART A: Abricate + VFDB ===" | tee -a "$LOG"

conda run -n clinical_env abricate \
    --db vfdb \
    --threads "$THREADS" \
    --minid 70 \
    --mincov 60 \
    "$ASSEMBLY_DIR"/*.fna \
    > "$OUT/abricate_vfdb_raw.tsv" 2>>"$LOG"

conda run -n clinical_env abricate \
    --summary "$OUT/abricate_vfdb_raw.tsv" \
    > "$OUT/abricate_vfdb_summary.tsv" 2>>"$LOG"

echo "[INFO] VFDB summary → $OUT/abricate_vfdb_summary.tsv" | tee -a "$LOG"

# ─── PART B: BLASTN vs local VFDB (mimics VFanalyzer) ────────────────────────
echo "[INFO] === PART B: BLAST vs local VFDB ===" | tee -a "$LOG"

VFDB_FASTA="$VFDB/VFDB_setB_nt.fas"
[[ ! -f "$VFDB_FASTA" ]] && VFDB_FASTA="$VFDB/vfdb.fna"
[[ ! -f "$VFDB_FASTA" ]] && VFDB_FASTA=$(find "$VFDB" -name "*.fna" -o -name "*.fas" 2>/dev/null | head -1 || true)

if [[ -n "${VFDB_FASTA:-}" && -f "$VFDB_FASTA" ]]; then
    VFDB_BLASTDB="$OUT/vfdb_blastdb"
    mkdir -p "$VFDB_BLASTDB"

    if [[ ! -f "$VFDB_BLASTDB/vfdb.nsi" ]]; then
        echo "[INFO] Building VFDB BLAST database..." | tee -a "$LOG"
        conda run -n clinical_env makeblastdb \
            -in "$VFDB_FASTA" -dbtype nucl \
            -out "$VFDB_BLASTDB/vfdb" -parse_seqids >> "$LOG" 2>&1
    fi

    BLAST_ALL="$OUT/blast_vfdb_all.tsv"
    echo -e "sample\tqseqid\tsseqid\tpident\tqcovs\tlength\tevalue\tbitscore\tgene_name" \
        > "$BLAST_ALL"

    while IFS= read -r SAMPLE; do
        GENOME="$ASSEMBLY_DIR/${SAMPLE}.fna"
        [[ ! -f "$GENOME" ]] && continue

        echo "[RUN ] BLAST: $SAMPLE vs VFDB" | tee -a "$LOG"

        conda run -n clinical_env blastn \
            -query "$GENOME" \
            -db "$VFDB_BLASTDB/vfdb" \
            -outfmt "6 qseqid sseqid pident qcovs length evalue bitscore stitle" \
            -perc_identity 70 -qcov_hsp_perc 60 \
            -num_threads "$THREADS" -evalue 1e-10 \
            2>>"$LOG" | \
        awk -F'\t' -v s="$SAMPLE" 'BEGIN{OFS="\t"} {
            gene = $2
            if (match($8, /\[[A-Za-z0-9_-]+\]/)) {
                gene = substr($8, RSTART+1, RLENGTH-2)
            }
            print s, $1, $2, $3, $4, $5, $6, $7, gene
        }' >> "$BLAST_ALL"
    done < "$SAMPLE_LIST"

    echo "[INFO] BLAST results → $BLAST_ALL" | tee -a "$LOG"
else
    echo "[WARN] VFDB FASTA not found at: ${VFDB_FASTA:-<empty>} — skipping BLAST step" | tee -a "$LOG"
fi

# ─── VFG CATEGORY MAP (A. baumannii specific) ────────────────────────────────
cat > "$OUT/vfg_category_map.tsv" << 'CATEOF'
gene	category
ompA	adherence
pilE	adherence
fliP	adherence
adeFGH	biofilm_formation
adeABC	biofilm_formation
csuA	biofilm_formation
csuB	biofilm_formation
csuC	biofilm_formation
csuD	biofilm_formation
csuE	biofilm_formation
pgaA	biofilm_formation
pgaB	biofilm_formation
pgaC	biofilm_formation
pgaD	biofilm_formation
bap	biofilm_formation
bfmR	biofilm_formation
bfmS	biofilm_formation
plc	enzyme
pld	enzyme
lpxA	immune_evasion
lpxB	immune_evasion
lpxC	immune_evasion
lpxD	immune_evasion
rmlD	immune_evasion
barA	iron_uptake
barB	iron_uptake
basA	iron_uptake
basB	iron_uptake
basC	iron_uptake
basD	iron_uptake
bauA	iron_uptake
bauB	iron_uptake
bauC	iron_uptake
bauD	iron_uptake
bauE	iron_uptake
bauF	iron_uptake
entE	iron_uptake
hemO	iron_uptake
abaI	regulation
abaR	regulation
katA	stress_adaptation
ompW	serum_resistance
omp33-36	serum_resistance
CATEOF

# ─── PRESENCE/ABSENCE MATRIX ─────────────────────────────────────────────────
PA_MATRIX="$RESULTS_DIR/13_summary_tables/VFG_presence_absence.tsv"
cp "$OUT/abricate_vfdb_summary.tsv" "$PA_MATRIX"

# ─── COUNT TABLE ─────────────────────────────────────────────────────────────
VFG_COUNT="$RESULTS_DIR/13_summary_tables/VFG_counts.tsv"
echo -e "sample\tn_VFG_total" > "$VFG_COUNT"

while IFS= read -r SAMPLE; do
    N=0
    if [[ -f "$OUT/abricate_vfdb_raw.tsv" ]]; then
        N=$(grep "/${SAMPLE}\.fna" "$OUT/abricate_vfdb_raw.tsv" \
            | awk -F'\t' '$10 >= 70 && $9 >= 60' | wc -l || echo 0)
    fi
    echo -e "${SAMPLE}\t${N}" >> "$VFG_COUNT"
done < "$SAMPLE_LIST"

echo "[INFO] VFG count table → $VFG_COUNT"
echo "[DONE] Step 04 complete — $(date)" | tee -a "$LOG"

echo ""
echo "=== Top 20 VFGs across all samples ===" | tee -a "$LOG"
awk -F'\t' 'NR>1 && $10>=70 {print $6}' "$OUT/abricate_vfdb_raw.tsv" \
    | sort | uniq -c | sort -rn | head -20 | tee -a "$LOG" || true
