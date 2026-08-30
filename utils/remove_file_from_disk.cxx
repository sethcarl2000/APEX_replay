
#include <APEX/utils.h>
// stdlib
#include <stdexcept>
#include <filesystem>
#include <system_error>

namespace APEX
{
namespace utils
{

void remove_file_from_disk(const std::string& path_str, PathCheckStatus* status)
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
    return;   
  }
  // claude opus 5 showed me how to do this (seth 29 aug 26)

  const fs::path path{path_str}; 

  std::error_code error_code; 

  //try to remove the file 
  const bool is_removed = fs::remove(path, error_code); 

  //check if there's an error code
  if (error_code) {
    if (status == nullptr) {
      throw std::system_error(error_code, "in <"+func+">: system error thrown when attempting to remove file: '"+ path.string()+"'");
    } else {
      status->message = "in <"+func+">: system error thrown when attempting to remove file: '"+ path.string()+"', message: " + error_code.message(); 
    }
    return; 
  }

  //now, let's check to make sure the path is absent. 
  if (fs::exists(path)) { 
    if (status == nullptr) {
      throw std::runtime_error("in <"+func+">: despite being deleted, path still exists: '"+path.string()+"' ");
    } else {
      status->message = "in <"+func+">: despite being deleted, path still exists: '"+path.string()+"'"; 
    }
    return;
  }
  
  if (status) status->message = "success"; 
}

}
}
