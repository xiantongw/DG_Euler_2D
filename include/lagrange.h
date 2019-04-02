#ifndef LAGRANGE
#define LAGRANGE

#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/qvm/mat_operations.hpp>

#include "../include/InvertMatrix.h"
#include "../include/utils.h"
#include "../include/GetQuadraturePointsWeight2D.h"
#include "../include/TriMesh.h"
#include "../include/ResData.h"
#include "../include/geometry.h"


namespace lagrange
{

    namespace ublas = boost::numeric::ublas;

    ublas::matrix<double> TriangleLagrange2D(int p);

    ublas::vector<ublas::vector<double> > MapReferenceToPhysicalLinear(ublas::vector<ublas::vector<double> > vertex, int p);

    ublas::vector<double> MapPhysicalToReferenceLinear(ublas::vector<ublas::vector<double> > vertex, ublas::vector<double> point, int p);

    ublas::vector<ublas::vector<double> > MapReferenceToPhysical(TriMesh mesh,  int ielem, int p, double (*pBumpFunction)(double));

    ublas::vector<ublas::matrix<double> > ConstructMassMatrix(int p, TriMesh mesh, ResData resdata);
    ublas::vector<ublas::matrix<double> > CalcInvMassMatrix(ublas::vector<ublas::matrix<double> > M);

    ublas::vector<double> CalcBaseFunction(ublas::matrix<double> TriLagrangeCoeff, double xi, double eta);

    ublas::matrix<double> CalcBaseFunctionGradient(ublas::matrix<double> TriLagrangeCoeff, double xi, double eta);

} // namespace lagrange

#endif