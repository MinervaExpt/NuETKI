#!/bin/bash

# ---- configuration ----
TODAY=$(date +"%b_%d")
#TODAY="Feb_25"
BASEDIR="/exp/minerva/app/users/cpernas/MAT_AL9/"
OUTBASE="/pnfs/minerva/scratch/users/cpernas/${TODAY}"

# output file lists (start clean)
CUT_LIST="CutSummary_files.txt"
DATA_LIST="DataHists_${TODAY}.txt"
MC_LIST="MCHists_${TODAY}.txt"

> "$CUT_LIST"
> "$DATA_LIST"
> "$MC_LIST"

# ---- playlists ----
playlists=(
  minervame1A
  minervame1B
  minervame1C
  minervame1D
  minervame1E
  minervame1F
  minervame1G
  minervame1L
  minervame1M
  minervame1N
  minervame1O
  minervame1P
)

# ---- submit jobs + record outputs ----
for pl in "${playlists[@]}"; do
  letter=$(echo "$pl" | sed -E 's/minervame1([A-Z])_.*/\1/')
  outdir="${OUTBASE}/${letter}"

  # submit job
  python SubmitJobsToGrid.py \
    --stage eventLoop \
    --playlist "$pl" \
    --basedir "$BASEDIR" \
    --outdir "$outdir" \
    --config "${BASEDIR}NuE_TKI/config/analysis.yaml"

  # record expected outputs
  echo "${outdir}/CutSummary.txt"        >> "$CUT_LIST"
  echo "${outdir}/Data_${TODAY}.root"    >> "$DATA_LIST"
  echo "${outdir}/MC_${TODAY}.root"      >> "$MC_LIST"
done

#copy the lists of output files to the base directory, helps to keep track of them
cp "$CUT_LIST" "${OUTBASE}"
cp "$DATA_LIST" "${OUTBASE}"
cp "$MC_LIST" "${OUTBASE}"
