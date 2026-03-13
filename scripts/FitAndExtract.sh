variables=(
  E_lep
  E_avail
  E_nu
  Lepton_Pt
  Lepton_Pl
  Theta_lep
  Proton_p
  Proton_Pt
  Theta_p
  Proton_T
  DeltaPt
  DeltaPtX
  DeltaPtY
  DeltaPl
  P_n
  AlphaPt
  PhiPt
)

for var in "${variables[@]}"; do
    cd /exp/minerva/data/users/cpernas/NuE_TKI/Feb_25/merged
    mkdir "$var"
    cd "$var"
    mkdir "xsec"
    FitBackgrounds ../data_merged.root ../mc_merged.root 1000000000 "$var"
    cd xsec
    ExtractCrossSection 1 ../../data_merged.root ../scaled_mc.root "$var"
done

cd /exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/scripts
