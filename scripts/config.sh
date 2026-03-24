#!/usr/bin/env bash
# =============================================================================
# config.sh — USER CONFIGURATION
# Edit the paths below to match your system before running the pipeline.
# This file is sourced by step00_setup.sh, which then auto-generates
# env_config.sh for use by all downstream steps.
# =============================================================================

# ─── PROJECT PATHS ───────────────────────────────────────────────────────────
# Root directory for this project (results, scripts, metadata will live here)
export PROJECT_DIR="$HOME/acinetobacter_genomics"

# Directory containing genome assemblies (.fna files, one per isolate)
export ASSEMBLY_DIR="$PROJECT_DIR/assembly"

# Root directory for all bioinformatics databases (see README for setup)
export DB_BASE="$HOME/databases/bacteria_microbial_genomics"

# ─── REFERENCE GENOME (Step 08 — SNP Phylogeny) ──────────────────────────────
# Path to A. baumannii ATCC 19606 reference genome (GCF_000015425.1)
# Download: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000015425.1/
export REFERENCE_GENOME="$PROJECT_DIR/references/ATCC_19606_genomic.fna"

# ─── COMPUTE SETTINGS ────────────────────────────────────────────────────────
export THREADS=8

# ─── DATABASE SUB-PATHS (relative to DB_BASE — change only if needed) ────────
export CARD_DB="$DB_BASE/card_db"
export VFDB="$DB_BASE/vfdb"
export GENOMAD_DB="$DB_BASE/genomad_db/genomad_db"
export MOB_DB="$DB_BASE/mob_suite_db"
export ABARI_DB="$DB_BASE/abari"
export GTDB_DATA="$DB_BASE/GTDB"
export AMRFINDER_DB="$DB_BASE/amrfinderplus_db/latest"
