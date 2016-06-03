#ifndef __SecDec_include_guard_%(name)s_integrands
#define __SecDec_include_guard_%(name)s_integrands
%(sector_includes)s

namespace %(name)s
{
    const unsigned int number_of_sectors = %(number_of_sectors)i;
    const unsigned int number_of_regulators = %(number_of_regulators)i;
    const unsigned int number_of_real_parameters = %(number_of_real_parameters)i;
    const unsigned int number_of_complex_parameters = %(number_of_complex_parameters)i;
    const std::array<int, number_of_regulators> requested_order = {%(requested_orders)s};
    const std::array<%(sectors_type)s, number_of_sectors> sectors = {%(sectors_initializer)s};
};
#endif
