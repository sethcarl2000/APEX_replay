
#include <APEX/utils.h>
// stdlib
#include <stdexcept>
#include <filesystem>
#include <system_error>
#include <cstdint> 

namespace APEX
{
namespace utils
{

double get_file_size_KB(const std::string& path_str, PathCheckStatus* status)
{
  namespace fs = std::filesystem; 

  //first, check if it's a regular file
  auto result_is_reg_file = is_path_regular_file(path_str); 
    
  const std::string func{__func__}; 

  if (auto result_is_reg_file = is_path_regular_file(path_str); !result_is_reg_file) {
    
    const std::string msg = "in <"+func+">: problem with 'regular_file' check. message: " + result_is_reg_file.message; 
    
    if (status == nullptr) {
      throw std::runtime_error(msg);
    } else {
      status->message = msg;
    }
    return utils::kNaN;   
  }
  
  const fs::path path{path_str}; 
  
  std::error_code error_code; 
  
  //try to remove the file 
  std::uintmax_t KB = fs::file_size(path, error_code);
  
  //check if there's an error code
  if (error_code) {
    if (status == nullptr) {
      throw std::system_error(error_code, "in <"+func+">: system error thrown when attempting to get file size: '"+ path.string()+"'");
    } else {
      status->message = "in <"+func+">: system error thrown when attempting to get file size: '"+ path.string()+"', message: " + error_code.message(); 
    }
    return utils::kNaN; 
  }
  if (status) status->message = "success";
  
  return (double)KB; 
}

}
}
