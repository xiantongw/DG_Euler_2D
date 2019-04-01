#include <iostream>
#include <algorithm>
#include <complex>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/minmax.hpp>

#include "../include/TriMesh.h"
#include "../include/utils.h"
#include "../include/lagrange.h"
#include "../include/geometry.h"
#include "../include/ConstructCurveMesh.h"
#include "../include/GetQuadraturePointsWeight2D.h"
#include "../include/solver.h"
#include "../include/euler.h"
#include "../include/Param.h"
#include "../include/Collective.h"
#include "../include/InvertMatrix.h"

using namespace std;
using namespace utils;
using namespace lagrange;


int main(int argc, char *argv[])
{
    namespace ublas = boost::numeric::ublas;

    Param param;
    // Set up the param struct
    param = ReadParamIn(string(argv[1]));
    TriMesh mesh(param.mesh_file);
    // Testing Calculate Residaul
    int p = param.order;
    int Np = int((p + 1) * (p + 2) / 2);
    TriMesh curved_mesh = mesh;
    string boundary_name="bottom";
    ConstructCurveMesh(mesh, curved_mesh, geometry::BumpFunction, boundary_name, p + 1);
    ublas::vector<double> States (curved_mesh.num_element * Np * 4, 0.0);
    // Construct free stream state
    ublas::vector<double> u_free = euler::CalcFreeStreamState_2DEuler(param);
    for (int ielem = 0; ielem < curved_mesh.num_element; ielem++)
    {
        for (int ip = 0; ip < Np; ip++)
        {
            States(ielem * Np * 4 + ip * 4 + 0) = u_free(0);
            States(ielem * Np * 4 + ip * 4 + 1) = u_free(1);
            States(ielem * Np * 4 + ip * 4 + 2) = u_free(2);
            States(ielem * Np * 4 + ip * 4 + 3) = u_free(3);
        }
    }
    ResData resdata = solver::CalcResData(curved_mesh, p);
    ublas::vector<double> dt(curved_mesh.E.size());
    ublas::vector<ublas::matrix<double> > M = lagrange::ConstructMassMatrix(p, curved_mesh, resdata);
    ublas::vector<ublas::matrix<double> > invM = lagrange::CalcInvMassMatrix(M);
    int MAXITER = param.MAXITER;
    int converged = 0;
    ublas::vector<double> States_new (curved_mesh.num_element * Np * 4, 0.0);
    for (int niter = 0; niter < MAXITER; niter++)
    {
        // cout << niter << endl;
        double norm_residual = 0.0;
        States_new = solver::TimeMarching_TVDRK3(curved_mesh, param, resdata, States, invM, p, converged, norm_residual);
	    if (niter % 100 == 0)
	    {
            std::cout << "NITER: " << niter << "\t" << "Residual Norm_Inf: ";
            cout.setf(ios::scientific, ios::floatfield);
            std::cout << setprecision(10) << norm_residual << std::endl;
        }
        States = States_new;
        if (converged)
           break;
    }
    return 0;
}
