#ifndef APEX_VDC_Hit_h
#define APEX_VDC_Hit_h

#include <APEX/EventHandler.h> 
#include <APEX/VDC.h> 

namespace APEX 
{
namespace VDC
{

/////////////////////////////////////////////////////////////////////////////////
class Hit { 

  bool   fis_RHRS; 
  int    fPlane;
  int    fWireNum; 
  double fRawTime; 
  double fT_offset; 

 public: 

  Hit( bool is_RHRS, int plane, int wire_num, double rawtime, double t_offset=0. ) 
    : fis_RHRS{is_RHRS}, 
    fPlane{plane}, 
    fWireNum{wire_num}, 
    fRawTime{rawtime},
    fT_offset{t_offset} {};  
  

  // methods which this class stores, i.e., this class' private data is the 'source of truth' 
  bool Is_RHRS() const { return fis_RHRS; } 
  int    Plane() const { return fPlane; } 
  int    WireNum() const { return fWireNum; } 
  double Rawtime() const { return fRawTime; }
  double T_offset() const { return fT_offset; }

  //set the timing offset for this hit. this is the only private data member we're allowed to modify once this class is already created. 
  void Set_T_offset(double t_offset) { fT_offset=t_offset; }

  // methods which this class fetches from tables, using its private raw data: 
  double Time() const { return VDC::RealTime(fis_RHRS,fPlane,fWireNum,fRawTime) - fT_offset; }
  double WirePos() const { return VDC::WirePos(fis_RHRS,fPlane,fWireNum); } 
  double W()     const { return VDC::w(fis_RHRS,fPlane); } 
  
}; 

}
}

#endif