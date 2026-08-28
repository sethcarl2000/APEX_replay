
#include <APEX/utils.h>

namespace APEX 
{
namespace utils 
{

template<> double intpow<2>(double x) { return x*x; }
template<> double intpow<3>(double x) { return x*x*x; }
template<> double intpow<4>(double x) { return x*x*x*x; }

}
}