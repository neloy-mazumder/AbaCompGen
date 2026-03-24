#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — MASTER PIPELINE
# A. baumannii Comparative Genomics
# Replication + Extension of Kumkar et al. 2022
#
# USAGE:
#   bash run_pipeline.sh              # run all steps
#   bash run_pipeline.sh --step 3     # run only step 3
#   bash run_pipeline.sh --from 5     # run steps 5 onwards
#   bash run_pipeline.sh --dry-run    # show steps without running
#
# PREREQUISITES:
#   1. Edit config.sh with your paths
#   2. Run step00_setup.sh at least once before using --step or --from
#
# PIPELINE OVERVIEW:
# ┌──────┬──────────────────────────────┬────────────────────────┐
# │ STEP │ SCRIPT                       │ CONDA ENV              │
# ├──────┼──────────────────────────────┼────────────────────────┤
# │  00  │ step00_setup.sh              │ base                   │
# │  01  │ step01_annotation.sh         │ annotation_env         │
# │  02  │ step02_mlst.sh               │ clinical_env           │
# │  03  │ step03_ARG_prediction.sh     │ clinical_env           │
# │  04  │ step04_VFG_prediction.sh     │ clinical_env           │
# │  05  │ step05_plasmid_genomad.sh    │ clinical_env           │
# │  06  │ step06_IS_isescan.sh         │ clinical_env           │
# │  07  │ step07_prophage.sh           │ clinical_env           │
# │  08  │ step08_SNP_phylogeny.sh      │ snippy2026 + phylo_env │
# │  09  │ step09_AbaRI_comM.sh         │ clinical_env           │
# │  10  │ step10_pangenome_roary.sh    │ roary_env              │ ← NOVEL
# │  11  │ step11_gtdbtk_temporal.sh    │ assembly_binning_env   │ ← NOVEL
# │  12  │ step12_summary_visualization │ r_env                  │
# └──────┴──────────────────────────────┴────────────────────────┘
#
# APPROXIMATE RUNTIMES (8 threads, ~47 genomes):
#   Step 01 (Prokka)      : ~20–40 min
#   Step 02 (MLST)        : ~2 min
#   Step 03 (ARG)         : ~10 min
#   Step 04 (VFG)         : ~15 min
#   Step 05 (Plasmids)    : ~45–90 min
#   Step 06 (IS)          : ~30 min
#   Step 07 (Prophage)    : ~5 min
#   Step 08 (Phylogeny)   : ~60–120 min
#   Step 09 (AbaRI)       : ~10 min
#   Step 10 (Pan-genome)  : ~20–40 min
#   Step 11 (GTDB)        : ~60–120 min
#   Step 12 (Summary)     : ~5 min
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─── LOAD CONFIG ─────────────────────────────────────────────────────────────
if [[ ! -f "$SCRIPT_DIR/config.sh" ]]; then
    echo "[ERROR] config.sh not found. Please copy config.sh.example → config.sh and set paths."
    exit 1
fi
source "$SCRIPT_DIR/config.sh"

# ─── PARSE ARGS ──────────────────────────────────────────────────────────────
ONLY_STEP=""
FROM_STEP=0
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --step)    ONLY_STEP="$2"; shift 2 ;;
        --from)    FROM_STEP="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ─── STEP REGISTRY ───────────────────────────────────────────────────────────
declare -a STEP_SCRIPTS=(
    "step00_setup.sh"
    "step01_annotation.sh"
    "step02_mlst.sh"
    "step03_ARG_prediction.sh"
    "step04_VFG_prediction.sh"
    "step05_plasmid_genomad.sh"
    "step06_IS_isescan.sh"
    "step07_prophage.sh"
    "step08_SNP_phylogeny.sh"
    "step09_AbaRI_comM.sh"
    "step10_pangenome_roary.sh"
    "step11_gtdbtk_temporal.sh"
    "step12_summary_visualization.sh"
)

# ─── RUN FUNCTION ────────────────────────────────────────────────────────────
run_step() {
    local IDX="$1"
    local SCRIPT="$2"
    local FULL_PATH="$SCRIPT_DIR/$SCRIPT"

    if [[ ! -f "$FULL_PATH" ]]; then
        echo "[ERROR] Script not found: $FULL_PATH"
        return 1
    fi

    printf "\n%s\n" "$(printf '=%.0s' {1..60})"
    echo "  RUNNING STEP $IDX: $SCRIPT"
    printf "%s\n\n" "$(printf '=%.0s' {1..60})"

    if $DRY_RUN; then
        echo "[DRY-RUN] Would execute: bash $FULL_PATH"
        return 0
    fi

    START_T=$(date +%s)
    bash "$FULL_PATH"
    END_T=$(date +%s)
    ELAPSED=$((END_T - START_T))
    printf "[TIMING] Step %02d completed in %d min %d sec\n" \
        "$IDX" "$((ELAPSED/60))" "$((ELAPSED%60))"
}

# ─── MAIN EXECUTION ──────────────────────────────────────────────────────────
PIPELINE_START=$(date +%s)

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   A. baumannii Comparative Genomics Pipeline             ║"
echo "║   Kumkar et al. 2022 — Replication + Extension           ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
echo "Script directory : $SCRIPT_DIR"
echo "Project dir      : $PROJECT_DIR"
echo "Mode             : $(${DRY_RUN} && echo 'DRY-RUN' || echo 'EXECUTE')"
[[ -n "$ONLY_STEP" ]] && echo "Single step      : $ONLY_STEP"
[[ "$FROM_STEP" -gt 0 ]] && echo "Starting from    : Step $FROM_STEP"
echo ""

for IDX in "${!STEP_SCRIPTS[@]}"; do
    SCRIPT="${STEP_SCRIPTS[$IDX]}"
    if [[ -n "$ONLY_STEP" ]]; then
        [[ "$IDX" != "$ONLY_STEP" ]] && continue
    elif [[ "$FROM_STEP" -gt 0 ]]; then
        [[ "$IDX" -lt "$FROM_STEP" ]] && continue
    fi
    run_step "$IDX" "$SCRIPT"
done

PIPELINE_END=$(date +%s)
TOTAL=$((PIPELINE_END - PIPELINE_START))

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   PIPELINE COMPLETE                                      ║"
printf "║   Total runtime: %d min %d sec%*s║\n" \
    "$((TOTAL/60))" "$((TOTAL%60))" "$((28-${#TOTAL}))" ""
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
echo "Key result directories:"
echo "  Annotation   : $RESULTS_DIR/01_annotation/"
echo "  MLST         : $RESULTS_DIR/02_mlst/"
echo "  ARGs         : $RESULTS_DIR/03_ARG_RGI/  +  $RESULTS_DIR/04_ARG_amrfinder/"
echo "  VFGs         : $RESULTS_DIR/05_VFG_abricate/"
echo "  Plasmids     : $RESULTS_DIR/06_plasmid_genomad/"
echo "  IS Elements  : $RESULTS_DIR/07_IS_isescan/"
echo "  Prophages    : $RESULTS_DIR/08_prophage_genomad/"
echo "  Phylogeny    : $RESULTS_DIR/09_SNP_phylogeny/iqtree/"
echo "  AbaRI        : $RESULTS_DIR/15_AbaRI_comM/"
echo "  Pan-genome   : $RESULTS_DIR/10_pangenome_roary/"
echo "  GTDB         : $RESULTS_DIR/12_gtdbtk/"
echo "  Master table : $RESULTS_DIR/13_summary_tables/master_genome_table.tsv"
echo "  Figures      : $RESULTS_DIR/figures/"
