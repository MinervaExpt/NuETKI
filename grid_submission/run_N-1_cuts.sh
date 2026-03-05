#!/bin/bash

# ---- configuration ----
TODAY=$(date +"%b_%d")
BASEDIR="/exp/minerva/app/users/cpernas/MAT_AL9/"
OUTBASE="/pnfs/minerva/scratch/users/cpernas/N-1_${TODAY}/"

# ---- N-1 plot config files ----
cut_plots=(
  StartPointVertexMultiplicity.yaml
  Afterpulsing.yaml
  DSCalVisE.yaml
  ODCalVisE.yaml
  VertexTrackMultiplicity.yaml
  TransverseGapScore.yaml
  NonMIPClusFrac.yaml
  EMLikeTrackScore.yaml
  MichelCut.yaml
  MeanFrontdEdX.yaml
  E_lep.yaml
  ModifiedEavailable.yaml
  NIsoBlobs.yaml
  ESC.yaml
  ExtraCuts.yaml
)

# ---- submit jobs + record outputs ----
for cut_config in "${cut_plots[@]}"; do
  cut_name="${cut_config%.yaml}"
  outdir="${OUTBASE}/${cut_name}"

  # submit job
  python SubmitJobsToGrid.py \
    --stage eventLoop \
    --playlist minervame1M \
    --basedir "$BASEDIR" \
    --outdir "$outdir" \
    --config "${BASEDIR}NuE_TKI/config/N-1/${cut_config}"

done
