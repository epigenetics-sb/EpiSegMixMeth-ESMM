#!/bin/bash

# Usage: ./merge_tracks.sh plus.ini minus.ini > merged.ini

PLUS_FILE="$1"
MINUS_FILE="$2"

if [[ -z "$PLUS_FILE" || -z "$MINUS_FILE" ]]; then
    echo "Usage: $0 <plus.ini> <minus.ini>" >&2
    exit 1
fi

# Extract file paths from each ini, preserving order
mapfile -t PLUS_FILES < <(grep '^file=' "$PLUS_FILE" | sed 's/file=//')
mapfile -t MINUS_FILES < <(grep '^file=' "$MINUS_FILE" | sed 's/file=//')

if [[ ${#PLUS_FILES[@]} -ne ${#MINUS_FILES[@]} ]]; then
    echo "Error: plus and minus files have different number of tracks (${#PLUS_FILES[@]} vs ${#MINUS_FILES[@]})" >&2
    exit 1
fi

# Extract the common track settings (everything except file=) from the plus file
# Assumes all tracks share the same settings
TRACK_SETTINGS=$(awk '/^\[RNA_Seq\]/{found=1; next} found && /^file=/{found=0; next} found{print}' "$PLUS_FILE" | sort -u)

for i in "${!PLUS_FILES[@]}"; do
    plus="${PLUS_FILES[$i]}"
    minus="${MINUS_FILES[$i]}"

    # Write plus strand track
    echo "[RNA_Seq]"
    echo "file=${plus}"
    echo "$TRACK_SETTINGS"
    echo ""

    # Write minus strand track (inverted)
    echo "[RNA_Seq]"
    echo "file=${minus}"
    echo "$TRACK_SETTINGS"
    echo "orientation = inverted"
    echo ""
done

# Append the spacer and gene track from the plus file (or minus — they're the same)
awk '/^\[spacer\]/,0' "$PLUS_FILE"
