#include <iostream>
#include <secdecutil/series.hpp>

int main()
{
    secdecutil::Series<int> exact(-2,1,{1,2,3,4},false,"eps");
    secdecutil::Series<int> truncated(-2,1,{1,2,3,4},true,"eps");
    secdecutil::Series<secdecutil::Series<int>> multivariate(1,2,
                                                             {
                                                                 {-2,-1,{1,2},false,"alpha"},
                                                                 {-2,-1,{3,4},false,"alpha"},
                                                             },false,"eps"
                                                             );

    std::cout << "exact:        " << exact << std::endl;
    std::cout << "truncated:    " << truncated << std::endl;
    std::cout << "multivariate: " << multivariate << std::endl << std::endl;

    std::cout << "exact + 1:         " << exact + 1 << std::endl;
    std::cout << "exact * exact:     " << exact * exact << std::endl;
    std::cout << "exact * truncated: " << exact * truncated << std::endl;
    std::cout << "exact.at(-2):      " << exact.at(-2) << std::endl;
}
