
#include <APEX/S2Event.h>
#include <APEX/units.h>
#include <APEX/run_parameters.h>
// stdlib
#include <stdexcept> 
#include <string> 

namespace APEX
{

//_________________________________________________________________________________________________________
S2Event::S2Event(bool is_RHRS, S2Hit hit) 
    : fis_RHRS{is_RHRS}, fS2Hit{hit} 
{
    if (!fS2Hit.IsCoinc()) {
        throw std::invalid_argument(
            "in <"+std::string{__func__}+">: S2 hit is not valid!"
        ); 
    }
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
void S2Event::make_vdc_hit_groups(
    int plane, 
    const ROOT::VecOps::RVec<double>& wire, 
    const ROOT::VecOps::RVec<double>& rawtime )
{
    if (plane < 0 || plane > VDC::kN_planes) {
        throw std::invalid_argument(
            "in <"+std::string{__func__}+">: VDC plane number '"+std::to_string(plane)+"' is not valid."
        ); 
        return; 
    }   

    auto& hit_groups = fHitGroups[plane]; 
    if (!hit_groups.empty()) {
        throw std::invalid_argument(
            "in <"+std::string{__func__}+">: VDC plane number '"+std::to_string(plane)+"' already has "+std::to_string(hit_groups.size())+" hit-groups made"
        ); 
        return; 
    }

    //make a buffer of all valid hits
    const auto nhits_total = wire.size(); 
    std::vector<VDC::Hit> buffer; buffer.reserve(nhits_total); 

    // the time of the S2 hit for this event.
    double t_S2 = S2Time(); 

    // make buffer of all valid VDC hits
    for (std::size_t i=0; i<nhits_total; i++) {

        VDC::Hit hit(Is_RHRS(), plane, static_cast<int>(wire[i]), rawtime[i], t_S2); 

        if (auto t = hit.Time(); 
            t > run_parameters::kVDC_min_realTime && 
            t < run_parameters::kVDC_max_realTime) 
        {
            buffer.emplace_back( std::move(hit) ); 
        }
    }

    //quit if there's not enough hits to form at least 1 hit. 
    if (buffer.empty() < run_parameters::kGroup_min_hits) return; 

    //now, we can make groups. 
    int prev_wire =-1; 

    size_t ind_begin=0; 
    size_t ind_end=ind_begin; 

    //add a new group by copying events from the current buffer 
    auto add_group = [&](){ 
        hit_groups.emplace_back(std::span<VDC::Hit>{
            buffer.begin()+ind_begin, 
            buffer.begin()+ind_end
        }); 
    };

    //keep iterating until we've gone through all hits
    for (size_t i=0; i<buffer.size(); i++) {

        const auto& hit = buffer[i]; 

#ifdef DEBUG_GROUP 
        printf(
            "plane %i\n"
            "   wire:       %i\n"
            "   time:       %.1f ns\n"
            "   raw time:   %.f\n",
            plane, 
            hit.WireNum(), 
            hit.Time()/units::ns, 
            hit.RawTime()
        );
#endif        

        int gap = ind_end==ind_begin ? 0 : hit.WireNum() - prev_wire -1; 
        if (gap > run_parameters::kGroup_max_gap) { // if the gap is too big to be part of one group 

            //does this cluster already have enough hits? 
            size_t n_hits_in_group = ind_end - ind_begin + 1; 
            if (n_hits_in_group >= run_parameters::kGroup_min_hits) {

                add_group(); 
            }     
            
            //otherwise, make a new group anyway.
            ind_begin = i;  

        }//if (span/gap too big) 

        //get ready to re-start this loop with the next hit
        prev_wire = hit.WireNum(); 
        ind_end = i; 
    }

    //We've reached the end, but we still might have enough hits to make a new group.  
    size_t n_hits_in_group = ind_end - ind_begin + 1; 
    if (n_hits_in_group >= run_parameters::kGroup_min_hits) {

        add_group(); 
    }     
}
//_________________________________________________________________________________________________________
//_________________________________________________________________________________________________________
//_________________________________________________________________________________________________________
//_________________________________________________________________________________________________________
//_________________________________________________________________________________________________________

}
