#ifndef APEX_PidManager_H
#define APEX_PidManager_H

#include <memory> 
#include <TObject.h> 
#include <APEX/PidDetector.h> 

/// @brief class which manages pointers to ApexPidDetector classes (helper classes used to process PID data)

namespace APEX
{

class PidManager {
public: 

    PidManager(); 

    //return the pre-shower detecotr (for LHRS: prl1)
    const PidDetector* GetPreShower(bool is_RHRS) const { return is_RHRS ? fDet_RHRS_ps.get() : fDet_LHRS_prl1.get(); } 

    //return the shower detector (for LHRS: prl2) 
    const PidDetector* GetShower(bool is_RHRS)    const { return is_RHRS ? fDet_RHRS_sh.get() : fDet_LHRS_prl2.get(); } 

private: 
    //left-arm shower detectorss
    std::unique_ptr<PidDetector> fDet_LHRS_prl1, fDet_LHRS_prl2; 

    //right-arm shower detectors
    std::unique_ptr<PidDetector> fDet_RHRS_sh, fDet_RHRS_ps; 
};

}

#endif