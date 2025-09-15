#!/bin/bash

# --- Configuration ---
SOURCE_DIR="/home/tfmuser/Desktop/3D_sims_alphabeta_rand/historical/data_offset_maxang=20"
DEST_DIR="/home/tfmuser/Desktop/3D_sims_alphabeta_rand/data"

# --- Script ---


echo "Copying files with index < 10..."


# Copy inc files
echo "  - Processing inc files..."
for file in "$SOURCE_DIR"/inc-*.npy; do
    number=$(basename "$file" .npy | cut -d'-' -f2)
    if (( number >= 10 && number < 30 )); then
        cp "$file" "$DEST_DIR"/
    fi
done




# Copy wavelength files
echo "  - Processing wavelength files..."
for file in "$SOURCE_DIR"/wvls-*.npy; do
    number=$(basename "$file" .npy | cut -d'-' -f2)
    if (( number >= 10 && number < 30 )); then
        cp "$file" "$DEST_DIR"/
    fi
done

# Copy incident angle files
echo "  - Processing incident angle files..."
for file in "$SOURCE_DIR"/inc_angles-*.npy; do
    number=$(basename "$file" .npy | cut -d'-' -f2)
    if (( number >= 10 && number < 30 )); then
        cp "$file" "$DEST_DIR"/
    fi
done

# Copy reflection data files
echo "  - Processing reflection files..."
for file in "$SOURCE_DIR"/refl-*.npy; do
    number=$(basename "$file" .npy | cut -d'-' -f2)
    if (( number >= 10 && number < 30 )); then
        cp "$file" "$DEST_DIR"/
    fi
done

echo "Done. Copied files are in the '$DEST_DIR' directory."
