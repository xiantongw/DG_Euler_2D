#ifndef CALCRESIDUAL_H
#define CALCRESIDUAL_H

#include <iostream>
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

// typedef struct ResData{
// 	int n_quad_1d;
// 	int n_quad_2d;
// 	int Np;
// 	int Nq;
// 	ublas::vector<double> x_quad_1d;
// 	ublas::vector<double> w_quad_1d;
// 	ublas::vector<double> x_quad_2d;
// 	ublas::vector<double> w_quad_2d;
// 	arr_2d Phi;
// 	arr_3d GPhi;
// 	arr_3d GPhi_Curved;
// 	arr_3d Phi_1D;
// 	arr_3d Phi_1D_Curved;
// 	arr_4d GPhi_1D;
// 	arr_4d GPhi_1D_Curved;
// } ResData;

ublas::vector<double> CalcResidual(TriMesh mesh, Param& param, ResData& resdata, ublas::vector<double> States, int p);

ResData CalcResData(TriMesh mesh, int p);

#endif