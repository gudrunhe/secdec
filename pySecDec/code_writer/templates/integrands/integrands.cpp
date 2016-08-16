#include <secdecutil/sector_container.hpp>
#include <secdecutil/series.hpp>
#include <vector>

#include "%(name)s.hpp"
%(sector_includes)s

namespace %(name)s
{
    const std::vector<%(sector_container_type)s> sectors = {%(sectors_initializer)s};
};
