#!/bin/bash

data_dir="/Users/shihab/Documents/FMRI/Control"
output_dir="/Users/shihab/Documents/NiLearn/Processed_FMRI/Control"
mkdir -p "$output_dir"

TR=3
mni_ref="$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz"

for func_in in "$data_dir"/*.nii.gz; do
    [[ ! -f "$func_in" ]] && continue

    echo "🔄 Processing: $func_in"
    base=$(basename "$func_in" .nii.gz)
    mc="${output_dir}/${base}_mc.nii.gz"
    st="${output_dir}/${base}_st.nii.gz"
    mean="${output_dir}/${base}_mean.nii.gz"
    mean_mni="${output_dir}/${base}_mean_mni.nii.gz"
    mat="${output_dir}/${base}_2mni.mat"
    smooth="${output_dir}/${base}_smooth.nii.gz"
    tmp_split_dir="${output_dir}/${base}_split"
    mkdir -p "$tmp_split_dir"

    # Step 1: Motion correction
    echo "📌 [1/6] Motion correction..."
    mcflirt -in "$func_in" -out "$mc" -refvol middle -plots || continue

    # Step 2: Slice timing correction
    echo "📌 [2/6] Slice timing correction..."
    slicetimer -i "$mc" -o "$st" --odd -r "$TR" || continue

    # Step 3: Compute mean image
    echo "📌 [3/6] Mean functional image..."
    fslmaths "$st" -Tmean "$mean" || continue

    # Step 4: Register mean to MNI
    echo "📌 [4/6] MNI registration..."
    flirt -in "$mean" -ref "$mni_ref" -out "$mean_mni" -omat "$mat" -dof 12 -cost corratio || continue

    # Step 5: Split and apply flirt to each volume
    echo "📌 [5/6] Apply transform to 4D (per volume)..."
    fslsplit "$st" "$tmp_split_dir/vol_" -t
    vol_list=()
    for vol in "$tmp_split_dir"/vol_*.nii.gz; do
        vol_out="${vol%.nii.gz}_mni.nii.gz"
        flirt -in "$vol" -ref "$mni_ref" -applyxfm -init "$mat" -out "$vol_out"
        vol_list+=("$vol_out")
    done

    # Merge volumes back into 4D
    merged="${output_dir}/${base}_mni.nii.gz"
    fslmerge -t "$merged" "${vol_list[@]}" || continue

    # Step 6: Smooth final 4D image
    echo "📌 [6/6] Smoothing..."
    fslmaths "$merged" -s 2 "$smooth" || continue

    echo "✅ Output: $smooth"

    # Optional: cleanup intermediate split volumes
    rm -rf "$tmp_split_dir"
done

