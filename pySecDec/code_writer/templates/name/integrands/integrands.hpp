#ifndef %(name)s_integrands_integrands_hpp_included
#define %(name)s_integrands_integrands_hpp_included
%(sector_includes)s

namespace %(name)s
{
    const unsigned int number_of_sectors = %(number_of_sectors)i;
    const unsigned int number_of_regulators = %(number_of_regulators)i;
    const unsigned int number_of_real_parameters = %(number_of_real_parameters)i;
    const unsigned int number_of_complex_parameters = %(number_of_complex_parameters)i;
    const std::array<int, number_of_regulators> lowest_orders = {%(lowest_orders)s}; // not including the prefactor
    const std::array<int, number_of_regulators> highest_orders = {%(highest_orders)s}; // not including the prefactor
    const std::array<int, number_of_regulators> requested_orders = {%(requested_orders)s};
    const std::array<%(integrand_container_type)s, number_of_sectors> sectors = {%(sectors_initializer)s};
};
#endif
