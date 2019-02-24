#ifndef UTILS
#define UTILS

#include <iostream>
#include <vector>

namespace utils{
    
    int getFullOrderIndex(int r, int s, int order);

    std::vector<int> getVertexIndex(std::vector<int>& element);

    bool sortcol_0(std::vector<int> const& v1, std::vector<int> const& v2);

    bool sortcol_2(std::vector<int> const& v1, std::vector<int> const& v2);

    template<typename T>
    std::vector<std::vector<T> > slice_by_row(std::vector<std::vector<T> > const& v, int m, int n);

    int missing_from_012(int input1, int input2);
    
}

#endif