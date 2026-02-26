#!/usr/bin/env bash
set -euo pipefail

# preprocess_fmri_fsl.sh
# ----------------------
# Minimal rs-fMRI preprocessing with FSL:
#   1) motion correction (mcflirt)
#   2) slice timing correction (slicetimer)
#   3) register mean image to MNI (flirt)
#   4) apply transform to each volume (split/apply/merge)
#   5) smoothing (fslmaths)
#
# This is a parameterized (path-flexible) version of the original script.
#
# Example:
#   bash scripts/fmri/preprocess_fmri_fsl.sh \
#     --input-dir /path/to/nifti_4d \
#     --output-dir /path/to/processed \
#     --tr 3 \
#     --smooth-sigma 2

usage() {
  cat <<EOF
Usage:
  preprocess_fmri_fsl.sh --input-dir DIR --output-dir DIR [options]

Required:
  --input-dir DIR       Folder with input 4D NIfTI files (*.nii.gz)
  --output-dir DIR      Folder to write processed outputs

Options:
  --tr FLOAT            Repetition time in seconds (default: 3)
  --mni-ref PATH        MNI reference brain (default: \$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz)
  --smooth-sigma FLOAT  Spatial smoothing sigma (default: 2)
  --slice-order STR     slicetimer flag: "odd" or "even" (default: odd)
  --pattern GLOB        File glob inside input-dir (default: *.nii.gz)
  --keep-split          Keep split volumes (default: delete)
  -h, --help            Show help

Notes:
  - Requires FSL commands: mcflirt, slicetimer, fslmaths, flirt, fslsplit, fslmerge
  - Output per input file: <base>_smooth.nii.gz
EOF
}

INPUT_DIR=""
OUTPUT_DIR=""
TR="3"
MNI_REF="${FSLDIR:-}/data/standard/MNI152_T1_2mm_brain.nii.gz"
SMOOTH_SIGMA="2"
SLICE_ORDER="odd"
PATTERN="*.nii.gz"
KEEP_SPLIT="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --tr) TR="$2"; shift 2 ;;
    --mni-ref) MNI_REF="$2"; shift 2 ;;
    --smooth-sigma) SMOOTH_SIGMA="$2"; shift 2 ;;
    --slice-order) SLICE_ORDER="$2"; shift 2 ;;
    --pattern) PATTERN="$2"; shift 2 ;;
    --keep-split) KEEP_SPLIT="1"; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1"; usage; exit 2 ;;
  esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
  echo "[ERROR] --input-dir and --output-dir are required."
  usage
  exit 2
fi

mkdir -p "$OUTPUT_DIR"

# Validate FSL tools
for cmd in mcflirt slicetimer fslmaths flirt fslsplit fslmerge; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "[ERROR] Missing command '$cmd'. Make sure FSL is installed and sourced."
    exit 3
  fi
done

if [[ ! -f "$MNI_REF" ]]; then
  echo "[ERROR] MNI reference not found: $MNI_REF"
  echo "Tip: pass --mni-ref explicitly, or ensure FSLDIR is set."
  exit 3
fi

shopt -s nullglob
files=("$INPUT_DIR"/$PATTERN)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "[ERROR] No files found in $INPUT_DIR with pattern '$PATTERN'"
  exit 4
fi

echo "[INFO] Input:  $INPUT_DIR"
echo "[INFO] Output: $OUTPUT_DIR"
echo "[INFO] TR=$TR | smooth_sigma=$SMOOTH_SIGMA | slice_order=$SLICE_ORDER"

for func_in in "${files[@]}"; do
  echo "🔄 Processing: $func_in"
  base=$(basename "$func_in" .nii.gz)

  mc="${OUTPUT_DIR}/${base}_mc.nii.gz"
  st="${OUTPUT_DIR}/${base}_st.nii.gz"
  mean="${OUTPUT_DIR}/${base}_mean.nii.gz"
  mean_mni="${OUTPUT_DIR}/${base}_mean_mni.nii.gz"
  mat="${OUTPUT_DIR}/${base}_2mni.mat"
  tmp_split_dir="${OUTPUT_DIR}/${base}_split"
  mkdir -p "$tmp_split_dir"

  # Step 1: Motion correction
  echo "📌 [1/6] Motion correction..."
  mcflirt -in "$func_in" -out "$mc" -refvol middle -plots

  # Step 2: Slice timing correction
  echo "📌 [2/6] Slice timing correction..."
  if [[ "$SLICE_ORDER" == "odd" ]]; then
    slicetimer -i "$mc" -o "$st" --odd -r "$TR"
  elif [[ "$SLICE_ORDER" == "even" ]]; then
    slicetimer -i "$mc" -o "$st" --even -r "$TR"
  else
    echo "[ERROR] --slice-order must be 'odd' or 'even' (got: $SLICE_ORDER)"
    exit 2
  fi

  # Step 3: Compute mean image
  echo "📌 [3/6] Mean functional image..."
  fslmaths "$st" -Tmean "$mean"

  # Step 4: Register mean to MNI
  echo "📌 [4/6] MNI registration..."
  flirt -in "$mean" -ref "$MNI_REF" -out "$mean_mni" -omat "$mat" -dof 12 -cost corratio

  # Step 5: Split and apply flirt to each volume
  echo "📌 [5/6] Apply transform to 4D (per volume)..."
  fslsplit "$st" "$tmp_split_dir/vol_" -t

  vol_list=()
  for vol in "$tmp_split_dir"/vol_*.nii.gz; do
    vol_out="${vol%.nii.gz}_mni.nii.gz"
    flirt -in "$vol" -ref "$MNI_REF" -applyxfm -init "$mat" -out "$vol_out"
    vol_list+=("$vol_out")
  done

  merged="${OUTPUT_DIR}/${base}_mni.nii.gz"
  fslmerge -t "$merged" "${vol_list[@]}"

  # Step 6: Smooth final 4D image
  echo "📌 [6/6] Smoothing..."
  smooth="${OUTPUT_DIR}/${base}_smooth.nii.gz"
  fslmaths "$merged" -s "$SMOOTH_SIGMA" "$smooth"

  echo "✅ Output: $smooth"

  if [[ "$KEEP_SPLIT" != "1" ]]; then
    rm -rf "$tmp_split_dir"
  fi
done
