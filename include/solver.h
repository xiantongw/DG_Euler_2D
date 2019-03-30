#ifndef CALCRESIDUAL_H
#define CALCRESIDUAL_H

#include <iostream>
#include <iomanip>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/multi_array.hpp>

#include "../include/InvertMatrix.h"
#include "../include/utils.h"
#include "../include/GetQuadraturePointsWeight1D.h"
#include "../include/GetQuadraturePointsWeight2D.h"
#include "../include/TriMesh.h"
#include "../include/lagrange.h"
#include "../include/geometry.h"
#include "../include/euler.h"
#include "../include/ConstructCurveMesh.h"
#include "../include/Param.h"
#include "../include/ResData.h"

namespace ublas = boost::numeric::ublas;

typedef boost::multi_array<double, 4> arr_4d;
typedef boost::multi_array<double, 3> arr_3d;
typedef boost::multi_array<double, 2> arr_2d;

namespace solver {

	ublas::vector<double> CalcResidual(TriMesh mesh, Param& param, ResData& resdata, ublas::vector<double> States, ublas::vector<double>& dtA, int p);

	ResData CalcResData(TriMesh mesh, int p);

	ublas::vector<double> TimeMarching(TriMesh mesh, Param& param, ResData& resdata, ublas::vector<double> States_old, ublas::vector<ublas::matrix<double> > invM, int p, int& converged);

}

#endif