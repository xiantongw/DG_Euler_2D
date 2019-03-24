#ifndef GEOMETRY
#define GEOMETRY

#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/qvm/mat_operations.hpp>

#include "../include/TriMesh.h"
#include "../include/utils.h"
#include "../include/lagrange.h"

namespace geometry
{
    namespace ublas = boost::numeric::ublas;

    double BumpFunction(double x);

    ublas::matrix<double> CalcJacobianLinear(TriMesh mesh, int ielem);
    ublas::matrix<double> CalcJacobianCurved(TriMesh mesh, int ielem, double xi, double eta);

    ublas::vector<int> GetEdgeLagrangeNodeIndex(int p, int local_index);
    ublas::vector<ublas::vector<double> > GetEdgeCoordinates(TriMesh mesh, int iedge);
    ublas::vector<int> GetEdgeCoordinatesIndex(TriMesh mesh, int iedge);
} // namespace geometry

#endif