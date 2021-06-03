#include "%(name)s.hpp"

#include <vector>
#include <memory> // std::shared_ptr, std::make_shared
#include <string>
#include <sstream>

#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation

#define INTEGRAL_NAME %(name)s

#include <secdecutil/pylink.hpp> // The python-C binding is general and therefore contained in the util
#include <secdecutil/pylink_amplitude.hpp>
