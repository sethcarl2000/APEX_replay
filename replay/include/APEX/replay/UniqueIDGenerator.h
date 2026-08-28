#ifndef APEX_replay_UniqueIDGenerator_h
#define APEX_replay_UniqueIDGenerator_h

#include <APEX/replay/EventCounter.h>

namespace APEX
{
namespace replay
{

class UniqueIDGenerator { 
  
 public: 
  /*** 
   *    Generates Unique ID each time it's called in thread-safe manner
   * 
   ***/

  UniqueIDGenerator() : fCounter{0} {};
  virtual ~UniqueIDGenerator() {}; 
  
  /// @return ID number which is unique across all threads 
  Long64_t GetID() { return fCounter.fetch_add(1, std::memory_order_relaxed); }

  //delete copy / move constructors
  UniqueIDGenerator(const UniqueIDGenerator&) = delete;
  UniqueIDGenerator& operator=(const UniqueIDGenerator&) = delete;
  
  
private:

  std::atomic<EventCounter> fCounter; 
  
}; 

}
}

#endif
