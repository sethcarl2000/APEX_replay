#ifndef APEX_utils_AtomicCounter_h
#define APEX_utils_AtomicCounter_h 

#include <atomic>
#include <memory> 
#include <string> 
#include <map> 
#include <mutex> 

class AtomicCounter {
private: 
     
    AtomicCounter(const AtomicCounter&) = delete; 
    AtomicCounter& operator=(const AtomicCounter&) = delete;  

    std::map<std::string, std::atomic<unsigned long long>> counters{}; 
    
    AtomicCounter() = default; 
    
public: 

    inline static AtomicCounter& Instance() {
        static AtomicCounter inst; 
        return inst; 
    }

    static void increment(std::string name) {

        // prevent races
        static std::mutex count_mutex; 
        std::lock_guard<std::mutex> g(count_mutex); 
        
        auto& counters = Instance().counters;
        
        auto find_it = counters.find(name); 
        if (find_it == counters.end()) {
            counters.emplace(name, 0); 
        } else {
            auto& count = find_it->second; 
            count.fetch_add(1, std::memory_order_relaxed);
        }
    }

    std::map<std::string, unsigned long long> static get() {
        auto& inst = Instance(); 
        if (inst.counters.empty()) return {}; 

        std::map<std::string, unsigned long long> counts;
        for(const auto& [name, counter] : inst.counters) {
            counts.emplace(name, static_cast<unsigned long long>(counter));
        }
        return counts; 
    }

}; 

#endif