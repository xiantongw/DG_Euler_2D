#ifndef EULER_H
#define EULER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include "../include/Param.h"
namespace euler
{
    namespace ublas = boost::numeric::ublas;

    ublas::matrix<double> CalcAnalyticalFlux(ublas::vector<double> state, double gamma);

    ublas::vector<double> CalcNumericalFlux(ublas::vector<double> uL, ublas::vector<double> uR, ublas::vector<double> norm,
                                            double gamma, char* type_flux, double& mws);

    ublas::vector<double> ApplyBoundaryCondition(ublas::vector<double> u, ublas::vector<double> norm,
                                                char* boundary_type, Param& cparam, double &mws);
}

#endif