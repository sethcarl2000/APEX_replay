/***
 * @class TapexVDCHitCluster.h
 * 
 * @brief Class which points to one VDC hit group, plus a computed intercept    
 * 
 * This 'hit cluster' marks a potential VDC intercept, and it also points to the TapexVDCHitGroup class
 * which all of it's hits belong to. 
 * 
 ***/

#ifndef TapexVDCHitCluster_H 
#define TapexVDCHitCluster_H

#include <TObject.h> 
#include <vector> 
#include "TapexVDCHitGroup.h"

class TapexVDCHitCluster : public TObject { 
  
 public: 
  TapexVDCHitCluster( TapexVDCHitGroup *group=nullptr, 
                      const double intercept=0, 
                      const double eta=0) 
    : fGroup{group}, fIntercept{intercept}, fEta_score{eta} {}; 
                  
  ~TapexVDCHitCluster() {/*noop*/}; 
  
  TapexVDCHitGroup* GetGroup() { return fGroup; }
  
  double Intercept() const { return fIntercept; }
  double Eta()       const { return fEta_score; }
  
 private: 
  TapexVDCHitGroup *fGroup; 
  double fIntercept; 
  double fEta_score; 
  
  ClassDef(TapexVDCHitCluster,0); 
};

#endif 