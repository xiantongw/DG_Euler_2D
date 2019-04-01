#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/qvm/mat_operations.hpp>

#include "./TriMesh.h"

namespace utils
{
    namespace ublas = boost::numeric::ublas;

    int GetFullOrderIndex(int r, int s, int order);

    std::vector<int> GetVertexIndex(std::vector<int>& element);

    bool SortByColumn0(std::vector<int> const& v1, std::vector<int> const& v2);

    bool SortByColumn2(std::vector<int> const& v1, std::vector<int> const& v2);

    template<typename T>
    std::vector<std::vector<T> > SliceByRow(std::vector<std::vector<T> > const& v, int m, int n);

    int MissingFrom012(int input1, int input2);

    template<typename T>
    ublas::vector<T> StdToBoostVector(std::vector<T> std_vec);

    double MaxBoostVector(ublas::vector<double> vec);

    ublas::matrix<double> Invert22Matrix(ublas::matrix<double> mat_input);

}

#endif