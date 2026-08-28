#ifndef APEX_run_parameters_H
#define APEX_run_parameters_H

#include <APEX/VDC.h> 

namespace APEX
{

/// @brief a basic namespace of constant parameters which effect the function of the VDC track reconstruction algorithm 
namespace run_parameters {

    constexpr int kVerbosity = 0; ///verbosity level. see functions/subroutines to see what different levels correspond to (default = 0)

    constexpr unsigned int kGroup_min_hits=2; ///minimum number of hits in a VDC group
    constexpr int kGroup_max_gap = 4;         ///maximum empty wire-gap in a VDC group

    constexpr double kVDC_max_realTime = 400e-9; ///maximum real-time for VDC hit, rel. to S2 hit
    constexpr double kVDC_min_realTime = -40e-9; ///minimum real-time for VDC hit, rel. to S2 hit

    /// this cut is how we determine if two adjacent S2 hits are 'twin' hits (likely caused by the same particle)
    constexpr double kCUT_S2_twinHit_timeErr =5e-9; 

    // spacing between coarse-grained grid-search 
    constexpr double kGridSpacing = ApexVDC::kWireSpacing/10.; 

    //min eta of one plane (during grid-searching). 
    // this is analogous to the num. of points found
    constexpr double kCUT_minEta = 1.950; 

    //track cuts
    constexpr double TRK_R_CUT_Dt     = 12e-9; 
    constexpr double TRK_L_CUT_Dt     = 20e-9; 
    constexpr double TRK_R_Dt_offset = +28e-9; 
    constexpr double TRK_L_Dt_offset = -15e-9;   
    constexpr double TRK_CUT_xParam = 0.60; 
    constexpr double TRK_R_xParam_offset = 0.10; 
    constexpr double TRK_L_xParam_offset = 0.00; 
    constexpr double TRK_CUT_Eta    = 6.; 
    constexpr double TRK_measureSigma = 5e-9; 
    constexpr int    TRK_CUT_nGoodPts_min_perPlane = 2; //good points per plane
    constexpr int    TRK_CUT_nGoodPts_min          = 10;   

    //for numerically-evaluated integrals with gaussians, this is the number of points sampled 
    constexpr int    kGausIntPoints = 20; 

    //these are kinematic cuts to the transport coordinate angles of VDC tracks. ph = dy/dz(tra), th = dx/dz(tra). 
    constexpr double CUT_ph_min[] = { -0.012 , -0.018 };
    constexpr double CUT_ph_max[] = {  0.012 ,  0.008 };   

    constexpr double CUT_th_min[] = { -0.018 , -0.020 };
    constexpr double CUT_th_max[] = {  0.018 ,  0.018 }; 
        
    constexpr double CUT_xParam = 2.00; 

    //window to cut on S2 coincidence hits 
    constexpr double CUT_twinHit_timeErr =5e-9; 
    
    // when there is a RHRS - LHRS S2 time conincidence, how many stddev's away from the center should we cut events? 
    constexpr double kS2_coinc_sigma_cut = 6.; 

    // if true, then multitrheadding is enabled 
    constexpr bool kEnableMT{true}; 

    // if 'true' then events with multiple coincidences will not be considered.
    constexpr bool kKillMultiCoincEvents{true}; 

    // minimum momentum to consider replaying a run (in MeV/c)
    constexpr double min_momentum = 800.; 

    // Range of S2 T_RHRS - T_LHRS to keep 
    constexpr double min_tr_tl = 10e-9; 
    constexpr double max_tr_tl = 60e-9; 
};

}

#endif
