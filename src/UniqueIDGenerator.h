#ifndef UniqueIDGenerator_H  
#define UniqueIDGenerator_H 

#include <atomic>
#include <EventCounter.h> 

class UniqueIDGenerator { 
  
 public: 
  /*** 
   *    Generates Unique ID each time it's called in thread-safe manner
   * 
   ***/

  UniqueIDGenerator() : fCounter{0} {};
  virtual ~UniqueIDGenerator() {}; 
  
  /// @return ID number which is unique across all threads 
  EventCounter GetID() { return fCounter.fetch_add(1, std::memory_order_relaxed); }

  //delete copy / move constructors
  UniqueIDGenerator(const UniqueIDGenerator&) = delete;
  UniqueIDGenerator& operator=(const UniqueIDGenerator&) = delete;
  
private:

  std::atomic<EventCounter> fCounter; 
  
}; 

#endif 
