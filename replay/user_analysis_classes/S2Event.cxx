
#include <APEX/S2Event.h>
#include <APEX/units.h>
#include <APEX/run_parameters.h>
// stdlib
#include <stdexcept> 
#include <string> 

namespace APEX
{

//_________________________________________________________________________________________________________
S2Event::S2Event(bool is_RHRS, S2Hit hit, const EventHandler* evt_handler, std::size_t _id) 
    : fis_RHRS{is_RHRS}, fS2Hit{hit}, id{_id}
{
    if (!fS2Hit.IsCoinc()) {
        throw std::invalid_argument(
            "in <"+std::string{__func__}+">: S2 hit is not valid!"
        ); 
    }

    if (!evt_handler) {
        throw std::invalid_argument(
            "in <"+std::string{__func__}+">: event handler passed is null"
        ); 
    }

    fEventHandler = std::make_shared<EventHandler>(
        Is_RHRS() ? replay::kRHRS : replay::kLHRS, 
        evt_handler->GetBeamCurrent(),
        evt_handler->GetRunNumber(), 
        Is_RHRS() ? &fS2Hit : evt_handler->GetS2Hit(false), 
        Is_RHRS() ? evt_handler->GetS2Hit(true)  : &fS2Hit   
    ); 
}
//_________________________________________________________________________________________________________
double S2Event::S2Time() const
{
    if (fS2Hit.Is_twinHit()) {
        // if it's a twin hit, then take the average of the two hits 
        return 0.5*( fS2Hit.Time() + fS2Hit.GetTwin()->Time() );  
    } else {
        return fS2Hit.Time(); 
    }
}
//_________________________________________________________________________________________________________
//_________________________________________________________________________________________________________
//_________________________________________________________________________________________________________
//_________________________________________________________________________________________________________
//_________________________________________________________________________________________________________
//_________________________________________________________________________________________________________

}
