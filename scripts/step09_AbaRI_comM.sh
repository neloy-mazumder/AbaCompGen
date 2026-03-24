#!/usr/bin/env bash
# =============================================================================
# STEP 09 — AbaR-type Resistance Island (AbaRI) Detection
# conda env : clinical_env  (blastn, makeblastdb, seqkit, entrez-direct)
#
# Method (matching Kumkar et al. 2022):
#   1. BLAST comM (A1S_1494 from ATCC 19606) against each genome
#   2. Extract ±30 kb flanking window around comM
#   3. BLAST window vs AbaRI island references (AbaR3, AbGRI1, AbaR4a)
#   4. BLAST window vs ARG marker sequences
#   5. Call AbaRI if island hit OR ≥2 ARG markers present
#   6. Classify variant from marker combination
#
# Note: Reference sequences are fetched from NCBI on first run.
#       Requires internet access on first run only.
# =============================================================================

set -euo pipefail
source "$PROJECT_DIR/env_config.sh"

OUT="$RESULTS_DIR/15_AbaRI_comM"
LOG="$RESULTS_DIR/00_logs/09_AbaRI.log"
FLANK_DIR="$OUT/flanking_regions"
BLAST_DIR="$OUT/blast_hits"
REF_DIR="$OUT/references"
mkdir -p "$FLANK_DIR" "$BLAST_DIR" "$REF_DIR" "$RESULTS_DIR/13_summary_tables"

WINDOW=30000   # ±30 kb flanking window

echo "[START] AbaRI detection — $(date)" | tee "$LOG"

# =============================================================================
# SECTION 1 — Build reference sequences (fetched from NCBI once)
# =============================================================================

# ── comM query (ATCC 19606 locus A1S_1494, ~1173 bp) ────────────────────────
comM_SEQ="$REF_DIR/comM_query.fna"

_comM_valid() {
    [[ -s "$comM_SEQ" ]] && \
    (( $(grep -v "^>" "$comM_SEQ" | tr -d '\n' | wc -c) >= 900 ))
}

if ! _comM_valid; then
    rm -f "$comM_SEQ"
    echo "[INFO] Fetching comM from NCBI (CP000521.1)..." | tee -a "$LOG"
    _TMP=$(mktemp /tmp/comM_cds_XXXXXX.fna)
    conda run -n clinical_env efetch \
        -db nuccore -id CP000521.1 -format fasta_cds_na > "$_TMP" 2>>"$LOG"
    awk 'found && /^>/{exit} /\[locus_tag=A1S_1494\]|\[gene=comM\]/{found=1} found{print}' \
        "$_TMP" > "$comM_SEQ"
    [[ ! -s "$comM_SEQ" ]] && \
        awk 'found && /^>/{exit} /gene=comM/{found=1} found{print}' \
        "$_TMP" > "$comM_SEQ"
    rm -f "$_TMP"
    if ! _comM_valid; then
        echo "[ERROR] comM fetch failed. Check network/entrez-direct install." | tee -a "$LOG"
        exit 1
    fi
    sed -i 's/^>.*/\>comM_ATCC19606_A1S_1494/' "$comM_SEQ"
    NT=$(grep -v "^>" "$comM_SEQ" | tr -d '\n' | wc -c)
    echo "[OK  ] comM query: ${NT} nt" | tee -a "$LOG"
fi

# ── AbaRI island reference sequences ─────────────────────────────────────────
ABARI_REF="$REF_DIR/AbaRI_island_refs.fna"
declare -A ABARI_ACC=( ["AbaR3"]="KF976462.1" ["AbGRI1"]="GQ900371.1" ["AbaR4a"]="KF976461.1" )

if [[ ! -s "$ABARI_REF" ]]; then
    > "$ABARI_REF"
    for VNAME in "${!ABARI_ACC[@]}"; do
        ACC="${ABARI_ACC[$VNAME]}"
        echo "[INFO] Fetching AbaRI reference: $VNAME ($ACC)" | tee -a "$LOG"
        _ATMP=$(mktemp /tmp/abari_XXXXXX.fna)
        conda run -n clinical_env efetch \
            -db nuccore -id "$ACC" -format fasta > "$_ATMP" 2>>"$LOG" || true
        if [[ -s "$_ATMP" ]]; then
            awk -v v="$VNAME" '/^>/{print ">"v"_reference"; next} {print}' \
                "$_ATMP" >> "$ABARI_REF"
            echo "[OK  ] $VNAME" | tee -a "$LOG"
        else
            echo "[WARN] Could not fetch $VNAME ($ACC)" | tee -a "$LOG"
        fi
        rm -f "$_ATMP"
    done
fi

# ── ARG marker sequences ──────────────────────────────────────────────────────
ARG_MARKERS="$REF_DIR/AbaRI_ARG_markers.fna"
declare -A ARG_ACC=(
    ["blaOXA-23"]="AY795964.1"  ["blaOXA-58"]="AY820284.1"
    ["tetB"]="AY534533.1"       ["strA"]="AF313472.1"
    ["sul2"]="AF324457.1"       ["aadA1"]="X12870.1"
    ["ISAba1"]="AJ277063.1"
)

if [[ ! -s "$ARG_MARKERS" ]]; then
    > "$ARG_MARKERS"
    for GENE in "${!ARG_ACC[@]}"; do
        ACC="${ARG_ACC[$GENE]}"
        _GTMP=$(mktemp /tmp/marker_XXXXXX.fna)
        conda run -n clinical_env efetch \
            -db nuccore -id "$ACC" -format fasta > "$_GTMP" 2>>"$LOG" || true
        if [[ -s "$_GTMP" ]]; then
            awk -v g="$GENE" '/^>/{print ">"g; next} {print}' "$_GTMP" >> "$ARG_MARKERS"
            echo "[OK  ] Marker: $GENE" | tee -a "$LOG"
        fi
        rm -f "$_GTMP"
    done
fi

for DB in "$ARG_MARKERS" "$ABARI_REF"; do
    if [[ -s "$DB" && ! -f "${DB}.nhr" ]]; then
        conda run -n clinical_env makeblastdb \
            -in "$DB" -dbtype nucl -out "$DB" >> "$LOG" 2>&1
    fi
done

# =============================================================================
# SECTION 2 — Per-sample analysis
# =============================================================================
ABARI_TABLE="$RESULTS_DIR/13_summary_tables/AbaRI_summary.tsv"
printf "sample\tcontig\tcomM_start\tcomM_end\tcomM_aln_len\tcomM_status\t"  > "$ABARI_TABLE"
printf "AbaRI_detected\tvariant\tisland_size_bp\t"                           >> "$ABARI_TABLE"
printf "blaOXA23\ttetB\tstrAB\tsul2\taadA\tISAba1\tn_markers\tn_ARGs\n"     >> "$ABARI_TABLE"

N_ABARI=0; N_TOTAL=0

while IFS= read -r SAMPLE; do
    GENOME="$ASSEMBLY_DIR/${SAMPLE}.fna"
    [[ ! -f "$GENOME" ]] && continue
    (( N_TOTAL++ ))
    echo "[RUN ] $SAMPLE" | tee -a "$LOG"

    COMHIT="$BLAST_DIR/${SAMPLE}_comM.tsv"
    conda run -n clinical_env blastn \
        -task blastn -query "$comM_SEQ" -subject "$GENOME" \
        -outfmt "6 sseqid sstart send length pident qcovs" \
        -perc_identity 75 -evalue 1e-5 \
        > "$COMHIT" 2>>"$LOG"

    if [[ ! -s "$COMHIT" ]]; then
        printf "%s\tNA\tNA\tNA\tNA\tNOT_FOUND\tno\tNA\tNA\t0\t0\t0\t0\t0\t0\t0\t0\n" \
            "$SAMPLE" >> "$ABARI_TABLE"
        echo "[INFO]  $SAMPLE — comM not found" | tee -a "$LOG"
        continue
    fi

    BEST=$(sort -t$'\t' -k4 -rn "$COMHIT" | head -1)
    CONTIG=$(echo "$BEST" | cut -f1)
    CS=$(echo "$BEST" | cut -f2); CE=$(echo "$BEST" | cut -f3)
    CLEN=$(echo "$BEST" | cut -f4)

    if (( CS > CE )); then S_MIN=$CE; S_MAX=$CS; else S_MIN=$CS; S_MAX=$CE; fi

    if   (( CLEN >= 1050 )); then COMSTATUS="intact"
    elif (( CLEN >=  500 )); then COMSTATUS="partial"
    else                          COMSTATUS="interrupted"
    fi

    CTGLEN=$(conda run -n clinical_env seqkit fx2tab -l -n "$GENOME" 2>/dev/null \
        | awk -F'\t' -v c="$CONTIG" '$1==c {print $2; exit}')
    [[ -z "$CTGLEN" || "$CTGLEN" -eq 0 ]] && CTGLEN=10000000

    FSTART=$(( S_MIN > WINDOW        ? S_MIN - WINDOW        : 1       ))
    FEND=$(( S_MAX + WINDOW < CTGLEN ? S_MAX + WINDOW        : CTGLEN  ))

    REGION="$FLANK_DIR/${SAMPLE}_flanking.fna"
    conda run -n clinical_env seqkit subseq \
        --chr "$CONTIG" -r "${FSTART}:${FEND}" \
        "$GENOME" > "$REGION" 2>>"$LOG" || cp "$GENOME" "$REGION"
    [[ ! -s "$REGION" ]] && cp "$GENOME" "$REGION"

    ISLAND_SIZE=0; VARIANT="none"
    if [[ -s "$ABARI_REF" && -f "${ABARI_REF}.nhr" ]]; then
        IBLAST="$BLAST_DIR/${SAMPLE}_island.tsv"
        conda run -n clinical_env blastn \
            -task blastn -query "$REGION" -db "$ABARI_REF" \
            -outfmt "6 sseqid pident length qcovs" \
            -perc_identity 80 -evalue 1e-20 \
            > "$IBLAST" 2>>"$LOG"
        if [[ -s "$IBLAST" ]]; then
            BEST_I=$(sort -t$'\t' -k3 -rn "$IBLAST" | head -1)
            ISLAND_SIZE=$(echo "$BEST_I" | cut -f3)
            VARIANT=$(echo "$BEST_I" | awk -F'\t' '{
                if ($1~/AbaR3/) print "AbaR3"
                else if ($1~/AbGRI/) print "AbGRI1"
                else if ($1~/AbaR4a/) print "AbaR4a"
                else print "novel"
            }')
        fi
    fi

    declare -A MH=( [blaOXA23]=0 [tetB]=0 [strAB]=0 [sul2]=0 [aadA]=0 [ISAba1]=0 )
    if [[ -s "$ARG_MARKERS" && -f "${ARG_MARKERS}.nhr" ]]; then
        MBLAST="$BLAST_DIR/${SAMPLE}_markers.tsv"
        conda run -n clinical_env blastn \
            -task blastn -query "$REGION" -db "$ARG_MARKERS" \
            -outfmt "6 sseqid pident length" \
            -perc_identity 80 -evalue 1e-10 \
            > "$MBLAST" 2>>"$LOG"
        if [[ -s "$MBLAST" ]]; then
            while IFS=$'\t' read -r MN PI ML; do
                case "$MN" in
                    blaOXA-23|blaOXA-58) MH[blaOXA23]=1 ;;
                    tetB)                MH[tetB]=1 ;;
                    strA|strB)           MH[strAB]=1 ;;
                    sul2)                MH[sul2]=1 ;;
                    aadA1)               MH[aadA]=1 ;;
                    ISAba1)              MH[ISAba1]=1 ;;
                esac
            done < "$MBLAST"
        fi
    fi

    N_M=$(( MH[blaOXA23] + MH[tetB] + MH[strAB] + MH[sul2] + MH[aadA] + MH[ISAba1] ))

    ABARI="no"
    if [[ "$VARIANT" != "none" ]] || (( N_M >= 2 )); then
        ABARI="yes"; (( N_ABARI++ ))
        if [[ "$VARIANT" == "none" ]]; then
            if   (( MH[blaOXA23]==1 && MH[tetB]==1 )); then VARIANT="AbaR3-like"
            elif (( MH[blaOXA23]==1 ));                 then VARIANT="AbaR4-like"
            elif (( MH[ISAba1]==1 ));                   then VARIANT="AbGRI-like"
            else                                             VARIANT="unclassified"
            fi
        fi
    fi

    N_ARG=0
    CARD_RAW="$RESULTS_DIR/03_ARG_RGI/abricate_card_raw.tsv"
    if [[ -f "$CARD_RAW" ]]; then
        N_ARG=$(awk -F'\t' -v samp="$SAMPLE" -v ctg="$CONTIG" \
            -v fs="$FSTART" -v fe="$FEND" '
            NR>1 {
                n=split($1,a,"/"); fname=a[n]; sub(/\.fna$/,"",fname)
                s=($4<$5)?$4:$5; e=($4>$5)?$4:$5
                if (fname==samp && $2==ctg && s>=fs && e<=fe) count++
            }
            END{print count+0}' "$CARD_RAW" 2>/dev/null || echo 0)
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SAMPLE" "$CONTIG" "$S_MIN" "$S_MAX" "$CLEN" "$COMSTATUS" \
        "$ABARI" "$VARIANT" "$ISLAND_SIZE" \
        "${MH[blaOXA23]}" "${MH[tetB]}" "${MH[strAB]}" \
        "${MH[sul2]}" "${MH[aadA]}" "${MH[ISAba1]}" \
        "$N_M" "$N_ARG" >> "$ABARI_TABLE"

    echo "[INFO]  $SAMPLE | comM: $COMSTATUS | AbaRI: $ABARI | variant: $VARIANT | markers: $N_M" \
        | tee -a "$LOG"
done < "$SAMPLE_LIST"

# =============================================================================
# SECTION 3 — Summary
# =============================================================================
echo "" | tee -a "$LOG"
echo "=== AbaRI Summary ===" | tee -a "$LOG"
printf "Genomes analyzed : %s\n" "$N_TOTAL"              | tee -a "$LOG"
printf "With AbaRI       : %s\n" "$N_ABARI"              | tee -a "$LOG"
printf "Without AbaRI    : %s\n" $(( N_TOTAL - N_ABARI)) | tee -a "$LOG"

echo "Variant breakdown:" | tee -a "$LOG"
awk -F'\t' 'NR>1 && $7=="yes"{v[$8]++}
    END{for(x in v) print "  "x": "v[x]}' "$ABARI_TABLE" | tee -a "$LOG"

echo "[INFO] Table → $ABARI_TABLE" | tee -a "$LOG"
echo "[DONE] AbaRI detection complete — $(date)" | tee -a "$LOG"
