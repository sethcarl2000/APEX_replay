#ifndef APEX_ObjRef_h
#define APEX_ObjRef_h 

#include <cstddef>
#include <limits> 
#include <vector> 
#include <map> 
#include <concepts>
#include <cstddef> 


namespace APEX
{

template<typename T> class ObjCollection; 

/// @brief This is a 'smart pointer', which let's analysis objects have reference to each other within the same S2 event 
/// @tparam T 
template<typename T> class ObjRef {
private: 

    // S2Event can handle the internal data from an ObjRef. 
    friend class ObjCollection<T>; 

    // this maximum possible value for an object's position will be used as the null value
    static constexpr std::size_t maxpos = std::numeric_limits<std::size_t>::max(); 

    // the index we want to point to. 
    std::size_t ind, event_id; 

    // only the S2Event can construct one of these refs that isn't null.  
    ObjRef(std::size_t _ind, std::size_t _event_id) : ind{_ind}, event_id{_event_id} {}; 

public: 

    /// @brief Anyone can use this default ctor (but the ref will be null.)
    ObjRef() : ind{maxpos}, event_id{maxpos} {}; 

    /// @return 'true' if this reference is null
    bool is_null() const { return ind==maxpos; }

    // this lets refs be keys of a std::map (sorted by ascending index number)
    bool operator<(const ObjRef& rhs) const { return ind < rhs.ind; }

};

}

#endif
