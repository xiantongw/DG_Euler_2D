#ifndef RESDATA_H
#define RESDATA_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/multi_array.hpp>

typedef boost::multi_array<double, 4> arr_4d;
typedef boost::multi_array<double, 3> arr_3d;
typedef boost::multi_array<double, 2> arr_2d;

typedef struct ResData{
	int n_quad_1d;
	int n_quad_2d;
	int Np;
	int Nq;
	ublas::vector<double> x_quad_1d;
	ublas::vector<double> w_quad_1d;
	ublas::vector<double> x_quad_2d;
	ublas::vector<double> w_quad_2d;
	arr_2d Phi;
	arr_3d GPhi;
	arr_3d GPhi_Curved;
	arr_3d Phi_1D;
	arr_3d Phi_1D_Curved;
	arr_4d GPhi_1D;
	arr_4d GPhi_1D_Curved;
	ublas::matrix<ublas::matrix<double> > jacobian_in_curved_elements;
	ublas::matrix<ublas::matrix<double> > invjacobian_in_curved_elements;
	ublas::vector<ublas::matrix<double> > jacobian_in_linear_elements;
	ublas::vector<ublas::matrix<double> > invjacobian_in_linear_elements;
} ResData;

#endif