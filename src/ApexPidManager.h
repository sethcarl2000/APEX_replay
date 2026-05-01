#ifndef ApexPidManager_H
#define ApexPidManager_H

#include <memory> 
#include <TObject.h> 
#include <ApexPidDetector.h> 

/// @brief class which manages pointers to ApexPidDetector classes (helper classes used to process PID data)

class ApexPidManager {
public: 

    ApexPidManager(); 

    virtual ~ApexPidManager() {}; 

    //return the pre-shower detecotr (for LHRS: prl1)
    const ApexPidDetector* GetPreShower(bool is_RHRS) const { return is_RHRS ? fDet_RHRS_ps.get() : fDet_LHRS_prl1.get(); } 

    //return the shower detector (for LHRS: prl2) 
    const ApexPidDetector* GetShower(bool is_RHRS)    const { return is_RHRS ? fDet_RHRS_sh.get() : fDet_LHRS_prl2.get(); } 

private: 
    //left-arm shower detectorss
    std::unique_ptr<ApexPidDetector> fDet_LHRS_prl1, fDet_LHRS_prl2; 

    //right-arm shower detectors
    std::unique_ptr<ApexPidDetector> fDet_RHRS_sh, fDet_RHRS_ps; 

    ClassDef(ApexPidManager,1); 
};

#endif