variables=(
  E_lep
  E_avail
  E_nu
  Lepton_Pt
  Lepton_Pl
  Theta_lep
  Proton_p
  #Proton_Pt #this one isn't filled in my GENIE things, an accident on my part
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
    cd "/exp/minerva/data/users/cpernas/NuE_TKI/Feb_25/merged/xsecs/closure_tests"
    root -b -q '/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/scripts/doClosureTest.C("../../'"${var}"'/xsec/'"${var}"'_crossSection.root", "../../GENIEXSECEXTRACT_all_playlists.root", "'"${var}"'")'
done

cd /exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/scripts
