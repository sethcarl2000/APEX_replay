#ifndef run_parameters_H
#define run_parameters_H


/// @brief a basic namespace of constant parameters which effect the function of the VDC track reconstruction algorithm 
namespace run_parameters {

    constexpr int kVerbosity = 0; ///verbosity level. see functions/subroutines to see what different levels correspond to (default = 0)

    constexpr unsigned int kGroup_min_hits=2; ///minimum number of hits in a VDC group
    constexpr int kGroup_max_gap = 4;         ///maximum empty wire-gap in a VDC group

    constexpr double kVDC_max_realTime = 400e-9; ///maximum real-time for VDC hit, rel. to S2 hit
    constexpr double kVDC_min_realTime = -40e-9; ///minimum real-time for VDC hit, rel. to S2 hit

    /// this cut is how we determine if two adjacent S2 hits are 'twin' hits (likely caused by the same particle)
    constexpr double kCUT_S2_twinHit_timeErr =5e-9; 

};

#endif