#ifndef replay_paths_h
#define replay_paths_h

#include <vector> 
#include <string>

namespace replay_paths
{
  constexpr size_t min_run = 3800;
  constexpr size_t max_run = 5000; 

  struct replay_segment {
    int min_run, max_run; 
    ULong64_t n_events{0};
    bool is_good{false};
    std::string path{"null"};
  }; 

  const std::vector<replay_segment> segment_laptop {
    {4175, 4199,  2313442,  true, "/home/sethcarl2000/j_research_linux/APEX_replay/replay-4175-4199.root"}
  }; 

  const std::vector<replay_segment> segments {
    {3800, 3824,        0, false, "/volatile/halla/apex/full_replay/production/replay-3800-3824.root"},
    {3825, 3849,        0, false, "/volatile/halla/apex/full_replay/production/replay-3825-3849.root"},
    {3850, 3874,     2500,  true, "/volatile/halla/apex/full_replay/production/replay-3850-3874.root"},
    {3875, 3899,   129801,  true, "/volatile/halla/apex/full_replay/production/replay-3875-3899.root"},
    {3900, 3924,    76328,  true, "/volatile/halla/apex/full_replay/production/replay-3900-3924.root"},
    {3925, 3949,   625131,  true, "/volatile/halla/apex/full_replay/production/replay-3925-3949.root"},
    {3950, 3974,   133925,  true, "/volatile/halla/apex/full_replay/production/replay-3950-3974.root"},
    {3975, 3999,  8294748,  true, "/volatile/halla/apex/full_replay/production/replay-3975-3999.root"},
    {4000, 4024,   644161,  true, "/volatile/halla/apex/full_replay/production/replay-4000-4024.root"},
    {4025, 4049,  1225406,  true, "/volatile/halla/apex/full_replay/production/replay-4025-4049.root"},
    {4050, 4074,   612399,  true, "/volatile/halla/apex/full_replay/production/replay-4050-4074.root"},
    {4075, 4099,   580927,  true, "/volatile/halla/apex/full_replay/production/replay-4075-4099.root"},
    {4100, 4124,        5, false, "/volatile/halla/apex/full_replay/production/replay-4100-4124.root"},
    {4125, 4149,        0, false, "/volatile/halla/apex/full_replay/production/replay-4125-4149.root"},
    {4150, 4174,        0, false, "/volatile/halla/apex/full_replay/production/replay-4150-4174.root"},
    {4175, 4199,  2313442,  true, "/volatile/halla/apex/full_replay/production/replay-4175-4199.root"},
    {4200, 4224, 38198197,  true, "/volatile/halla/apex/full_replay/production/replay-4200-4224.root"},
    {4225, 4249, 37966339,  true, "/volatile/halla/apex/full_replay/production/replay-4225-4249.root"},
    {4250, 4274, 37318859,  true, "/volatile/halla/apex/full_replay/production/replay-4250-4274.root"},
    {4275, 4299, 26729053,  true, "/volatile/halla/apex/full_replay/production/replay-4275-4299.root"},
    {4300, 4324, 32925084,  true, "/volatile/halla/apex/full_replay/production/replay-4300-4324.root"},
    {4325, 4349, 30958517,  true, "/volatile/halla/apex/full_replay/production/replay-4325-4349.root"},
    {4350, 4374, 11728688,  true, "/volatile/halla/apex/full_replay/production/replay-4350-4374.root"},
    {4375, 4399, 37436364,  true, "/volatile/halla/apex/full_replay/production/replay-4375-4399.root"},
    {4400, 4424, 13804469,  true, "/volatile/halla/apex/full_replay/production/replay-4400-4424.root"},
    {4425, 4449, 36903691,  true, "/volatile/halla/apex/full_replay/production/replay-4425-4449.root"},
    {4450, 4474, 32585915,  true, "/volatile/halla/apex/full_replay/production/replay-4450-4474.root"},
    {4475, 4499,   530378,  true, "/volatile/halla/apex/full_replay/production/replay-4475-4499.root"},
    {4500, 4524,  1902208,  true, "/volatile/halla/apex/full_replay/production/replay-4500-4524.root"},
    {4525, 4549,  1038107,  true, "/volatile/halla/apex/full_replay/production/replay-4525-4549.root"},
    {4550, 4574,  1801767,  true, "/volatile/halla/apex/full_replay/production/replay-4550-4574.root"},
    {4575, 4599,   706348,  true, "/volatile/halla/apex/full_replay/production/replay-4575-4599.root"},
    {4600, 4624,  1399328,  true, "/volatile/halla/apex/full_replay/production/replay-4600-4624.root"},
    {4625, 4649,        0, false, "/volatile/halla/apex/full_replay/production/replay-4625-4649.root"},
    {4650, 4674,   621088,  true, "/volatile/halla/apex/full_replay/production/replay-4650-4674.root"},
    {4675, 4699,  1724694,  true, "/volatile/halla/apex/full_replay/production/replay-4675-4699.root"},
    {4700, 4724,  1515897,  true, "/volatile/halla/apex/full_replay/production/replay-4700-4724.root"},
    {4725, 4749,        0, false, "/volatile/halla/apex/full_replay/production/replay-4725-4749.root"},
    {4750, 4774,        0, false, "/volatile/halla/apex/full_replay/production/replay-4750-4774.root"},
    {4775, 4799,        0, false, "/volatile/halla/apex/full_replay/production/replay-4775-4799.root"},
    {4800, 4824,        0, false, "/volatile/halla/apex/full_replay/production/replay-4800-4824.root"},
    {4825, 4849,        0, false, "/volatile/halla/apex/full_replay/production/replay-4825-4849.root"},
    {4850, 4874,        0, false, "/volatile/halla/apex/full_replay/production/replay-4850-4874.root"},
    {4875, 4899,        0, false, "/volatile/halla/apex/full_replay/production/replay-4875-4899.root"},
    {4900, 4924,        0, false, "/volatile/halla/apex/full_replay/production/replay-4900-4924.root"},
    {4925, 4949,        0, false, "/volatile/halla/apex/full_replay/production/replay-4925-4949.root"},
    {4950, 4974,        0, false, "/volatile/halla/apex/full_replay/production/replay-4950-4974.root"},
    {4975, 4999,        0, false, "/volatile/halla/apex/full_replay/production/replay-4975-4999.root"}
  }; 
  
  //list of all replay files (most recent replay)
  const std::vector<std::string> list{
    "/volatile/halla/apex/full_replay/production/replay-3800-3824.root",
    "/volatile/halla/apex/full_replay/production/replay-3825-3849.root",
    "/volatile/halla/apex/full_replay/production/replay-3850-3874.root",
    "/volatile/halla/apex/full_replay/production/replay-3875-3899.root",
    "/volatile/halla/apex/full_replay/production/replay-3900-3924.root",
    "/volatile/halla/apex/full_replay/production/replay-3925-3949.root",
    "/volatile/halla/apex/full_replay/production/replay-3950-3974.root",
    "/volatile/halla/apex/full_replay/production/replay-3975-3999.root",
    "/volatile/halla/apex/full_replay/production/replay-4000-4024.root",
    "/volatile/halla/apex/full_replay/production/replay-4025-4049.root",
    "/volatile/halla/apex/full_replay/production/replay-4050-4074.root",
    "/volatile/halla/apex/full_replay/production/replay-4075-4099.root",
    "/volatile/halla/apex/full_replay/production/replay-4100-4124.root",
    "/volatile/halla/apex/full_replay/production/replay-4125-4149.root",
    "/volatile/halla/apex/full_replay/production/replay-4150-4174.root",
    "/volatile/halla/apex/full_replay/production/replay-4175-4199.root",
    "/volatile/halla/apex/full_replay/production/replay-4200-4224.root",
    "/volatile/halla/apex/full_replay/production/replay-4225-4249.root",
    "/volatile/halla/apex/full_replay/production/replay-4250-4274.root",
    "/volatile/halla/apex/full_replay/production/replay-4275-4299.root",
    "/volatile/halla/apex/full_replay/production/replay-4300-4324.root",
    "/volatile/halla/apex/full_replay/production/replay-4325-4349.root",
    "/volatile/halla/apex/full_replay/production/replay-4350-4374.root",
    "/volatile/halla/apex/full_replay/production/replay-4375-4399.root",
    "/volatile/halla/apex/full_replay/production/replay-4400-4424.root",
    "/volatile/halla/apex/full_replay/production/replay-4425-4449.root",
    "/volatile/halla/apex/full_replay/production/replay-4450-4474.root",
    "/volatile/halla/apex/full_replay/production/replay-4475-4499.root",
    "/volatile/halla/apex/full_replay/production/replay-4500-4524.root",
    "/volatile/halla/apex/full_replay/production/replay-4525-4549.root",
    "/volatile/halla/apex/full_replay/production/replay-4550-4574.root",
    "/volatile/halla/apex/full_replay/production/replay-4575-4599.root",
    "/volatile/halla/apex/full_replay/production/replay-4600-4624.root",
    "/volatile/halla/apex/full_replay/production/replay-4625-4649.root",
    "/volatile/halla/apex/full_replay/production/replay-4650-4674.root",
    "/volatile/halla/apex/full_replay/production/replay-4675-4699.root",
    "/volatile/halla/apex/full_replay/production/replay-4700-4724.root",
    "/volatile/halla/apex/full_replay/production/replay-4725-4749.root",
    "/volatile/halla/apex/full_replay/production/replay-4750-4774.root",
    "/volatile/halla/apex/full_replay/production/replay-4775-4799.root",
    "/volatile/halla/apex/full_replay/production/replay-4800-4824.root",
    "/volatile/halla/apex/full_replay/production/replay-4825-4849.root",
    "/volatile/halla/apex/full_replay/production/replay-4850-4874.root",
    "/volatile/halla/apex/full_replay/production/replay-4875-4899.root",
    "/volatile/halla/apex/full_replay/production/replay-4900-4924.root",
    "/volatile/halla/apex/full_replay/production/replay-4925-4949.root",
    "/volatile/halla/apex/full_replay/production/replay-4950-4974.root",
    "/volatile/halla/apex/full_replay/production/replay-4975-4999.root"//*/ 
  };
  
}; 

#endif
