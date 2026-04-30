#ifndef BINNING_H
#define BINNING_H

#include <vector>

//quick lambda function for scaling vector, i'm using to convert from GeV to MeV (or back)
inline std::vector<double> scale(std::vector<double> v, double factor) {
    for (auto& x : v) x *= factor;
    return v;
}
// ============================================================================
// Analysis Variables
// =============================================================================

inline const std::vector<double> electronEnergyBins = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,20};
inline const std::vector<double> EavailBins = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1, 1.2};
inline const std::vector<double> leptonPt_bins = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
inline const std::vector<double> Pt_bins = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
inline const std::vector<double> Pl_bins = {0, 3, 4.5, 5.5, 6.5, 8, 12, 30};
inline const std::vector<double> electronAngleBins = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40};

inline const std::vector<double> protonMomentumBins = {0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375};
inline const std::vector<double> protonPtBins = {0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.1};
inline const std::vector<double> protonAngleBins = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180};
inline const std::vector<double> T_p_bins = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; 

//inline const std::vector<double> deltaPt_bins = {0, 0.2, 0.4, 0.6, 1, 1.5, 3};
//These are the delta pt bins used in Jeffrey's published TKI paper
inline const std::vector<double> deltaPt_bins = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.4};
//inline const std::vector<double> deltaPtXBins = {-2, -1.8, -1.6, -1.4, -1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2};
inline const std::vector<double> deltaPtXBins = {-2.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4, 0.7, 2.0};
//Jeffrey's again
inline const std::vector<double> deltaPtYBins = {-2.0, -1.1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 1.3};
//Copied from Jeffrey
//inline const std::vector<double> deltaPlBins = {-1.0, -0.25, 0.00, 0.10, 0.20, 0.30, 0.40, 0.60, 1.0};
inline const std::vector<double> deltaPlBins = {-0.25, 0.00, 0.10, 0.20, 0.30, 0.40, 0.60};
inline const std::vector<double> PnBins = {0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 1.0};
//Jeffrey has bins twice this width, so like 0, 40, 80, etc. I don't think I need that?
inline const std::vector<double> alphaAngleBins = {0,20,40,60,80,100,120,140,160,180};
//inline const std::vector<double> phiAngleBins = {0,10,20,30,40,50,60,70,80,90,100};
//These are from Jeffrey's paper:
inline const std::vector<double> phiAngleBins = {0,10,30,50,80,120,180};


inline const std::vector<double> KE_bins = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0, 1.25, 1.5};

//pt mu (dan): [0, 0.15, 0.25, 0.33, 0.40, 0.47, 0.55, 0.7, 0.85, 1, 1.25, 2.50]
//pt z (jeffrey): {0, 3, 4.5, 5.5, 6.5, 8, 12, 30}
//pt mu (jeffrey): {0, 0.1, 0.2, 0.3, 0.45, 0.75, 1.25, 5} 
//pt e (Hang): {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6}
//delta pt (jeffrey): {0, 0.2, 0.4, 0.6, 1, 1.5, 2.5} 

// ============================================================================
// Cut quantities
// =============================================================================
inline const std::vector<double> binaryCut_bins = {-0.5,0.5,1.5}; //For anything that's either 0 or 1

//tracker region: minZ =5980, maxZ = 8422,
inline const std::vector<double> VertexZ_bins {5880, 5980, 6080, 6180, 6280, 6380, 6480, 6580, 6680, 6780, 6880, 6980, 7080, 7180, 7280, 7380, 7480, 7580, 7680, 7780, 7880, 7980, 8080, 8180, 8280, 8380, 8480};

inline const std::vector<double> StartPointVertexMultiplicity_bins = {-1,0,1,2,3,4,5,6,7,8};

inline const std::vector<double> Afterpulsing_bins = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1};

inline const std::vector<double> Deadtime_bins = {0,1,2,3,4,5,6};

inline const std::vector<double> DSCalVisE_bins = {0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25};

inline const std::vector<double> ODCalVisE_bins = {0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04,0.0425,0.045,0.0475,0.05,0.0525,0.055,0.0575,0.06};

inline const std::vector<double> VertexTrackMultiplicity_bins = {0,1,2,3,4,5,6};

//inline const std::vector<double> TransverseGapScore_bins = {0,2,4,6,8,10,12,14,15,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50};
inline const std::vector<double> TransverseGapScore_bins = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};

inline const std::vector<double> NonMIPClusFrac_bins = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1};

//inline const std::vector<double> EMScore_bins = {0.499,0.549,0.599,0.649,0.699,0.749,0.799,0.849,0.899,0.949,0.999,1.049};
inline const std::vector<double> EMScore_bins = {0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05};

inline const std::vector<double> NMichel_bins = {0, 1, 2, 3, 4, 5};

inline const std::vector<double> MeanFrontDEDX_bins = {0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,8,9,10};

//inline const std::vector<double> Modified_E_avail_bins = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2};
inline const std::vector<double> Modified_E_avail_bins = scale({0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2}, 1000);

inline const std::vector<double> ESCChi2_bins = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40};

inline const std::vector<double> Psi_bins = {0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};

inline const std::vector<double> Etheta_bins = {0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.015,0.02,0.025};

//inline const std::vector<double> E_lep_Sin2Theta_bins = {0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.015,0.02,0.025,0.03};
//inline const std::vector<double> Etheta_bins = {0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};

#endif
