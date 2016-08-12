#ifndef %(name)s_integrands_integrands_hpp_included
#define %(name)s_integrands_integrands_hpp_included
%(sector_includes)s

namespace %(name)s
{
    const unsigned int number_of_sectors = %(number_of_sectors)i;
    const unsigned int number_of_regulators = %(number_of_regulators)i; //TODO: names of regulators
    const unsigned int number_of_real_parameters = %(number_of_real_parameters)i; //TODO: names of real_parameters
    const unsigned int number_of_complex_parameters = %(number_of_complex_parameters)i; //TODO: names of complex_parameters
    const std::vector<int> lowest_orders = {%(lowest_orders)s}; // not including the prefactor // TODO: lowest_prefactor_orders
    const std::vector<int> highest_orders = {%(highest_orders)s}; // not including the prefactor // TODO: highest_prefactor_orders
    const std::vector<int> requested_orders = {%(requested_orders)s};
    const std::vector<%(sector_container_type)s> sectors = {%(sectors_initializer)s};
    // TODO: prefactor
};
#endif
