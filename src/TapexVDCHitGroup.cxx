#include "TapexVDCHitGroup.h"
#include <stdexcept> 
#include <sstream> 
#include <limits> 

namespace {
    //a NAN double to return upon arrival at invalid funciton state
    constexpr double return_nan = std::numeric_limits<double>::quiet_NaN(); 
};      

//____________________________________________________________________________________________________________________
int TapexVDCHitGroup::WireNum(unsigned int h) const 
{
#ifdef RANGE_CHECKS
    if (h >= (unsigned int)fHits.size()) {
        std::ostringstream oss; 
        oss << "in <TapexVDCHitGroup::" << __func__ << ">: Invalid hit-number (" << h << ")."; 
        throw std::invalid_argument(oss.str()); 
        return -1; 
    }
#endif
    return fHits[h].wNum(); 
}
//____________________________________________________________________________________________________________________
double TapexVDCHitGroup::WirePos(unsigned int h) const 
{
#ifdef RANGE_CHECKS
    if (h >= (unsigned int)fHits.size()) {
        std::ostringstream oss; 
        oss << "in <TapexVDCHitGroup::" << __func__ << ">: Invalid hit-number (" << h << ")."; 
        throw std::invalid_argument(oss.str()); 
        return return_nan; 
    }
#endif
    return fHits[h].wPos(); 
}
//____________________________________________________________________________________________________________________
double TapexVDCHitGroup::Time(unsigned int h) const 
{
#ifdef RANGE_CHECKS
    if (h >= (unsigned int)fHits.size()) {
        std::ostringstream oss; 
        oss << "in <TapexVDCHitGroup::" << __func__ << ">: Invalid hit-number (" << h << ")."; 
        throw std::invalid_argument(oss.str()); 
        return return_nan; 
    }
#endif
    return fHits[h].Time(); 
}
//____________________________________________________________________________________________________________________
TapexVDCHit& TapexVDCHitGroup::GetHit(unsigned int h)
{
#ifdef RANGE_CHECKS
    if (h >= (unsigned int)fHits.size()) {
        std::ostringstream oss; 
        oss << "in <TapexVDCHitGroup::" << __func__ << ">: Invalid hit-number (" << h << ")."; 
        throw std::invalid_argument(oss.str()); 
        return fHits[0]; 
    }
#endif
    return fHits[h]; 
}
//____________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________

ClassImp(TapexVDCHitGroup); 