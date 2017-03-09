#include <iostream>
#include <secdecutil/series.hpp>

int main()
{
    secdecutil::Series<int> exact(-2,1,{1,2,3,4},false,"eps");
    secdecutil::Series<int> truncated(-2,1,{1,2,3,4},true,"eps");

    std::cout << "exact:     " << exact << std::endl;
    std::cout << "truncated: " << truncated << std::endl << std::endl;

    std::cout << "exact + 1:         " << exact + 1 << std::endl;
    std::cout << "exact * exact:     " << exact * exact << std::endl;
    std::cout << "exact * truncated: " << exact * truncated << std::endl;
    std::cout << "exact.at(-2):      " << exact.at(-2) << std::endl;
}
