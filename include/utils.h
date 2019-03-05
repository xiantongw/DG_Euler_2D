#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <cmath>

namespace utils
{
    int GetFullOrderIndex(int r, int s, int order);

    std::vector<int> GetVertexIndex(std::vector<int>& element);

    bool SortByColumn0(std::vector<int> const& v1, std::vector<int> const& v2);

    bool SortByColumn2(std::vector<int> const& v1, std::vector<int> const& v2);

    template<typename T>
    std::vector<std::vector<T> > SliceByRow(std::vector<std::vector<T> > const& v, int m, int n);

    int MissingFrom012(int input1, int input2);

}

#endif