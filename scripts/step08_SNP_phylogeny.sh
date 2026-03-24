#!/usr/bin/env bash
# =============================================================================
# STEP 08 — SNP-Based Phylogenetic Analysis
# conda env : snippy2026  (snippy 4.0.2, snp-sites, snp-dists)
#             phylogeny_env  (IQ-TREE 2.x, FastTree, MAFFT)
# Tools     : snippy       — pairwise SNP calling vs reference
#             snippy-core  — core SNP alignment
#             IQ-TREE 2    — maximum-likelihood phylogeny (1000 UFBoot)
#             FastTree     — rapid pre-visualization tree
# Paper ref : "CSIPhylogeny — A. baumannii ATCC 19606 as reference"
# Reference : A. baumannii ATCC 19606 (GCF_000015425.1)
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

PHYLO_DIR="$RESULTS_DIR/09_SNP_phylogeny"
SNIPPY_OUT="$PHYLO_DIR/snippy_out"
CORE_DIR="$PHYLO_DIR/core_snps"
IQTREE_DIR="$PHYLO_DIR/iqtree"
LOG="$RESULTS_DIR/00_logs/08_phylogeny.log"
mkdir -p "$SNIPPY_OUT" "$CORE_DIR" "$IQTREE_DIR"

echo "[INFO] Starting SNP Phylogeny — $(date)" | tee "$LOG"

# ─── VALIDATE REFERENCE GENOME ───────────────────────────────────────────────
if [[ ! -f "$REFERENCE_GENOME" ]]; then
    echo "[ERROR] Reference genome not found: $REFERENCE_GENOME" | tee -a "$LOG"
    echo "        Set REFERENCE_GENOME in config.sh and download:" | tee -a "$LOG"
    echo "        https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000015425.1/" | tee -a "$LOG"
    exit 1
fi
echo "[OK  ] Reference: $REFERENCE_GENOME" | tee -a "$LOG"

# ─── PATCH: snippy 4.0.2 samtools version check bug ─────────────────────────
# snippy uses Perl string comparison for version numbers, causing "1.23" < "1.7"
# (lexicographic). Fix: replace with proper numeric version check.
SNIPPY_BIN=$(conda run -n snippy2026 which snippy 2>/dev/null || true)
if [[ -n "$SNIPPY_BIN" && -f "$SNIPPY_BIN" ]]; then
    if grep -q 'version->parse.*ge.*version->parse' "$SNIPPY_BIN" 2>/dev/null; then
        echo "[FIX ] Patching snippy samtools version check in $SNIPPY_BIN" | tee -a "$LOG"
        perl -i -pe '
            s{version->parse\(\$ver\)\s+ge\s+version->parse\(\$MIN_SAMTOOLS\)}
             {do { my \@a=split(\/\\.\/,\$ver); my \@b=split(\/\\.\/,\$MIN_SAMTOOLS); \$a[0]*10000+(\$a[1]\/\/0) >= \$b[0]*10000+(\$b[1]\/\/0) }}g
        ' "$SNIPPY_BIN"
        echo "[OK  ] snippy version check patched." | tee -a "$LOG"
    fi
fi

# ─── RESOLVE IQ-TREE BINARY ──────────────────────────────────────────────────
if conda run -n phylogeny_env iqtree2 --version >> "$LOG" 2>&1; then
    IQTREE_BIN="iqtree2"
elif conda run -n phylogeny_env iqtree --version >> "$LOG" 2>&1; then
    IQTREE_BIN="iqtree"
else
    echo "[WARN] No iqtree binary found in phylogeny_env — tree step will be skipped." | tee -a "$LOG"
    IQTREE_BIN=""
fi
echo "[OK  ] IQ-TREE binary: ${IQTREE_BIN:-not found}" | tee -a "$LOG"

# ─── STEP 8.1: snippy per sample ─────────────────────────────────────────────
while IFS= read -r SAMPLE; do
    GENOME="$ASSEMBLY_DIR/${SAMPLE}.fna"
    [[ ! -f "$GENOME" ]] && continue

    if [[ -d "$SNIPPY_OUT/$SAMPLE" && -f "$SNIPPY_OUT/$SAMPLE/snps.vcf" ]]; then
        echo "[SKIP] snippy: $SAMPLE" | tee -a "$LOG"
        continue
    fi

    echo "[RUN ] snippy: $SAMPLE" | tee -a "$LOG"
    conda run -n snippy2026 snippy \
        --outdir "$SNIPPY_OUT/$SAMPLE" \
        --ref "$REFERENCE_GENOME" \
        --ctgs "$GENOME" \
        --cpus "$THREADS" \
        --force \
        >> "$LOG" 2>&1 || \
    echo "[WARN] snippy failed for $SAMPLE" | tee -a "$LOG"
done < "$SAMPLE_LIST"

# ─── STEP 8.2: snippy-core ───────────────────────────────────────────────────
echo "[INFO] Running snippy-core..." | tee -a "$LOG"

SNIPPY_DIR_ARGS=()
while IFS= read -r S; do
    if [[ -d "$SNIPPY_OUT/$S" && -f "$SNIPPY_OUT/$S/snps.vcf" ]]; then
        SNIPPY_DIR_ARGS+=("$SNIPPY_OUT/$S")
    fi
done < "$SAMPLE_LIST"

if [[ ${#SNIPPY_DIR_ARGS[@]} -eq 0 ]]; then
    echo "[ERROR] No completed snippy results — cannot run snippy-core." | tee -a "$LOG"
    exit 1
fi

cd "$CORE_DIR"
conda run -n snippy2026 snippy-core \
    --ref "$REFERENCE_GENOME" \
    --prefix core \
    "${SNIPPY_DIR_ARGS[@]}" \
    >> "$LOG" 2>&1

echo "[INFO] Core SNP alignment → $CORE_DIR/core.aln" | tee -a "$LOG"

# ─── STEP 8.3: Remove constant sites ─────────────────────────────────────────
echo "[INFO] Stripping constant sites with snp-sites..." | tee -a "$LOG"
conda run -n snippy2026 snp-sites \
    -c "$CORE_DIR/core.aln" -o "$CORE_DIR/core.snp-only.aln" \
    >> "$LOG" 2>&1 || \
conda run -n phylogeny_env snp-sites \
    -c "$CORE_DIR/core.aln" -o "$CORE_DIR/core.snp-only.aln" \
    >> "$LOG" 2>&1 || \
{ echo "[WARN] snp-sites not found — using full core alignment"
  cp "$CORE_DIR/core.aln" "$CORE_DIR/core.snp-only.aln"; }

# ─── STEP 8.4: IQ-TREE ML phylogeny ─────────────────────────────────────────
echo "[INFO] Running IQ-TREE (GTR+F+ASC, 1000 UFBoot)..." | tee -a "$LOG"
if [[ -n "$IQTREE_BIN" ]]; then
    conda run -n phylogeny_env "$IQTREE_BIN" \
        -s "$CORE_DIR/core.snp-only.aln" \
        -m GTR+F+ASC \
        -B 1000 \
        -T "$THREADS" \
        --prefix "$IQTREE_DIR/acinetobacter_snp" \
        --redo \
        >> "$LOG" 2>&1
    echo "[INFO] IQ-TREE done → $IQTREE_DIR/acinetobacter_snp.treefile" | tee -a "$LOG"
else
    echo "[WARN] Skipping IQ-TREE — binary not found." | tee -a "$LOG"
fi

# ─── STEP 8.5: FastTree (rapid pre-viz) ──────────────────────────────────────
conda run -n phylogeny_env FastTree \
    -nt -gtr "$CORE_DIR/core.snp-only.aln" \
    > "$IQTREE_DIR/acinetobacter_fasttree.nwk" 2>>"$LOG" || \
echo "[WARN] FastTree not available" | tee -a "$LOG"

# ─── STEP 8.6: iTOL annotation files ─────────────────────────────────────────
ITOL_DIR="$IQTREE_DIR/iTOL_annotations"
mkdir -p "$ITOL_DIR"

cat > "$ITOL_DIR/ST_colors.txt" << 'ITOLEOF'
DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	Sequence Type
COLOR	#ff0000
DATA
ITOLEOF

awk -F'\t' 'NR>1{
    if ($2 == "2")       col="#e41a1c"
    else if ($2 == "622") col="#377eb8"
    else if ($2 == "10")  col="#4daf4a"
    else if ($2 == "25")  col="#984ea3"
    else if ($2 == "1")   col="#ff7f00"
    else col="#999999"
    print $1"\t"col"\tST"$2
}' "$RESULTS_DIR/02_mlst/mlst_pasteur.tsv" \
    >> "$ITOL_DIR/ST_colors.txt" 2>/dev/null || true

cat > "$ITOL_DIR/metadata_labels.txt" << 'LABELEOF'
DATASET_TEXT
SEPARATOR TAB
DATASET_LABEL	Hospital_Year
COLOR	#000000
DATA
LABELEOF
awk -F'\t' 'NR>1{print $1"\t"$3"_"$2"\t-1\t#333333\tnormal\t1\t0"}' \
    "$METADATA" >> "$ITOL_DIR/metadata_labels.txt" 2>/dev/null || true

echo "[INFO] iTOL files → $ITOL_DIR"
echo "[INFO] Upload tree to https://itol.embl.de/ and drag annotation files"

# ─── STEP 8.7: SNP distance matrix ───────────────────────────────────────────
conda run -n snippy2026 snp-dists \
    "$CORE_DIR/core.snp-only.aln" \
    > "$CORE_DIR/snp_distance_matrix.tsv" 2>>"$LOG" || \
echo "[WARN] snp-dists not available — skipping distance matrix" | tee -a "$LOG"

echo "[DONE] Step 08 complete — $(date)" | tee -a "$LOG"
echo ""
echo "Key outputs:"
echo "  Tree (IQ-TREE) : $IQTREE_DIR/acinetobacter_snp.treefile"
echo "  SNP matrix     : $CORE_DIR/snp_distance_matrix.tsv"
echo "  iTOL annotation: $ITOL_DIR/"
