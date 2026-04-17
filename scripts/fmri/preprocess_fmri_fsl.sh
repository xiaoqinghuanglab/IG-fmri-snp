#!/bin/bash

# ======================================
# SLURM DIRECTIVES (Resource Allocation)
# ======================================
#SBATCH --job-name=ADNI3_Year5
#SBATCH -A r01042
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --array=1-99
#SBATCH --output=IG-fmri-snp/ADNI_Outputs/derivatives/logs/%x/fmriprep_%A_%a.out
#SBATCH --error=IG-fmri-snp/ADNI_Outputs/derivatives/logs/%x/fmriprep_%A_%a.err

# Load modules
module purge
module load apptainer

# =============================================================================
# CONFIGURATION & PATHS
# =============================================================================
DATASET_NAME="Data_folder_Name"
BASE_DIR="IG-fmri-snp"
BIDS_DIR="${BASE_DIR}/ADNI_Outputs/datasets/${DATASET_NAME}"
DERIVATIVES_DIR="${BASE_DIR}/ADNI_Outputs/derivatives/${DATASET_NAME}"
WORK_DIR="/N/scratch/${USER}/fmriprep_work/job_${SLURM_ARRAY_JOB_ID}_task_${SLURM_ARRAY_TASK_ID}"
TF_HOME="${BASE_DIR}/ADNI_Outputs/work/templateflow"
FS_LICENSE="${BASE_DIR}/Hasan/license.txt"
FMRIPREP_SIF="${BASE_DIR}/Hasan/fmriprep-25.2.3.sif"

FS_SUBJECTS_DIR="${DERIVATIVES_DIR}/sourcedata/freesurfer"

# =============================================================================
# FSAVERAGE TEMPLATE SETUP (Race-condition safe)
# Task 1 runs setup; all others wait until it's done.
# =============================================================================
TEMPLATE_READY_FLAG="${FS_SUBJECTS_DIR}/.templates_ready"

if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    echo "Task 1: Running fsaverage template setup..."

    mkdir -p "${FS_SUBJECTS_DIR}"

    # Clean stale/incomplete templates
    rm -rf "${FS_SUBJECTS_DIR}/fsaverage"
    rm -rf "${FS_SUBJECTS_DIR}/fsaverage5"
    rm -f  "${TEMPLATE_READY_FLAG}"

    apptainer exec --cleanenv \
        -B "${FS_SUBJECTS_DIR}":/out \
        "${FMRIPREP_SIF}" \
        cp -r /opt/freesurfer/subjects/fsaverage /out/

    apptainer exec --cleanenv \
        -B "${FS_SUBJECTS_DIR}":/out \
        "${FMRIPREP_SIF}" \
        cp -r /opt/freesurfer/subjects/fsaverage5 /out/

    # Signal to all other tasks that templates are ready
    touch "${TEMPLATE_READY_FLAG}"
    echo "Task 1: Template setup complete. Flag written."

else
    echo "Task ${SLURM_ARRAY_TASK_ID}: Waiting for template setup (Task 1)..."
    WAIT_SECONDS=0
    MAX_WAIT=300  # 5 minute timeout

    while [ ! -f "${TEMPLATE_READY_FLAG}" ]; do
        sleep 10
        WAIT_SECONDS=$((WAIT_SECONDS + 10))
        if [ "${WAIT_SECONDS}" -ge "${MAX_WAIT}" ]; then
            echo "ERROR: Timed out waiting for template flag after ${MAX_WAIT}s. Exiting."
            exit 1
        fi
    done

    echo "Task ${SLURM_ARRAY_TASK_ID}: Templates ready. Proceeding."
fi

# =============================================================================
# PIPELINE EXECUTION
# =============================================================================
if [ ! -f "${FS_LICENSE}" ]; then
    echo "ERROR: FreeSurfer license not found at ${FS_LICENSE}"
    exit 1
fi

SUBJECTS_FILE="${BIDS_DIR}/participants.tsv"
RAW_SUB=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SUBJECTS_FILE}" | awk '{print $1}')

if [ -z "${RAW_SUB}" ]; then
    echo "Subject not found on line $((SLURM_ARRAY_TASK_ID + 1)) of TSV. Exiting."
    exit 0
fi

BIDS_ID=${RAW_SUB#sub-}
FS_OUTPUT="${FS_SUBJECTS_DIR}/sub-${BIDS_ID}"

if [ -d "${FS_OUTPUT}" ]; then
    echo "==================================================="
    echo "SKIPPING: ${RAW_SUB}"
    echo "Reason: FreeSurfer output already exists at:"
    echo "  ${FS_OUTPUT}"
    echo "==================================================="
    exit 0
fi

BOLD_CHECK=$(find "${BIDS_DIR}/${RAW_SUB}/func" -name "*_bold.nii.gz" -print -quit 2>/dev/null)

if [ -z "${BOLD_CHECK}" ]; then
    echo "WARNING: No BOLD data found for ${RAW_SUB}. Switching to ANATOMICAL-ONLY mode."
    EXTRA_FLAGS="--anat-only"
else
    echo "SUCCESS: BOLD data found for ${RAW_SUB}. Running full pipeline."
    EXTRA_FLAGS=""
fi

mkdir -p "${DERIVATIVES_DIR}" "${WORK_DIR}" "${TF_HOME}"

echo "Processing: ${RAW_SUB}"
echo "Mode:       ${EXTRA_FLAGS:-Full Pipeline}"

# =============================================================================
# EXECUTION
# =============================================================================
apptainer run --cleanenv \
    -B ${BIDS_DIR}:/data:ro \
    -B ${DERIVATIVES_DIR}:/out \
    -B ${WORK_DIR}:/work \
    -B ${TF_HOME}:/templateflow \
    -B ${FS_LICENSE}:/opt/freesurfer/license.txt \
    --env TEMPLATEFLOW_HOME=/templateflow \
    ${FMRIPREP_SIF} \
    /data /out participant \
    --participant-label ${BIDS_ID} \
    --work-dir /work \
    --nthreads 8 \
    --omp-nthreads 4 \
    --mem-mb 90000 \
    --output-spaces MNI152NLin2009cAsym:res-2 anat func fsaverage5 \
    --output-layout bids \
    --use-syn-sdc warn \
    --cifti-output \
    --skull-strip-t1w auto \
    --bids-database-dir /work/bids_db \
    --skip-bids-validation \
    --notrack \
    ${EXTRA_FLAGS}

# ============================================================================
# POST-PROCESSING CLEANUP
# ============================================================================
if [ $? -eq 0 ]; then
    echo "SUCCESS: sub-${BIDS_ID} completed. Removing work directory."
    rm -rf "${WORK_DIR}"
else
    echo "FAILURE: sub-${BIDS_ID} failed. Check logs and work directory at ${WORK_DIR}"
    exit 1
fi
