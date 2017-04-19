#include <iostream>
#include <secdecutil/integrand_container.hpp>

int main()
{
    using input_t = const double * const;
    using return_t = double;

    std::function<return_t(input_t)> f1 = [] (input_t x) { return 2*x[0]; };
    secdecutil::IntegrandContainer<return_t,input_t> c1(1,f1);

    std::function<return_t(input_t)> f2 = [] (input_t x) { return x[0]*x[1]; };
    secdecutil::IntegrandContainer<return_t,input_t> c2(2,f2);

    secdecutil::IntegrandContainer<return_t,input_t> c3 = c1 + c2;
    const double point[]{1.0,2.0};

    std::cout << "c1.number_of_integration_variables: " << c1.number_of_integration_variables << std::endl;
    std::cout << "c2.number_of_integration_variables: " << c2.number_of_integration_variables << std::endl << std::endl;

    std::cout << "c3.number_of_integration_variables: " << c3.number_of_integration_variables << std::endl;
    std::cout << "c3.integrand(point):                " << c3.integrand(point) << std::endl;
}
