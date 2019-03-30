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
    for (int ielem = 0; ielem < curved_mesh.num_element; ielem++)
    {
        for (int ip = 0; ip < Np; ip++)
        {
            States(ielem * Np * 4 + ip * 4 + 0) = 1.0;
            States(ielem * Np * 4 + ip * 4 + 1) = param.mach_inf * cos(param.attack_angle);
            States(ielem * Np * 4 + ip * 4 + 2) = param.mach_inf * sin(param.attack_angle);
            States(ielem * Np * 4 + ip * 4 + 3) = 1 / (param.gamma * (param.gamma - 1)) +
                                                    0.5 * param.mach_inf * param.mach_inf;
        }
    }
    ResData resdata = solver::CalcResData(curved_mesh, p);
    ublas::vector<double> dtA(curved_mesh.E.size(), 0.0);
    ublas::vector<double> Residual = solver::CalcResidual(curved_mesh, param, resdata, States, dtA, p);
    cout << "Total Residual Norm: " << ublas::norm_2(Residual) << endl;
    ublas::vector<ublas::matrix<double> > M = lagrange::ConstructMassMatrix(p, curved_mesh, resdata);
    ublas::vector<ublas::matrix<double> > invM = lagrange::CalcInvMassMatrix(M);

    int MAXITER = 10000000;
    int converged = 1;
    for (int niter = 0; niter < MAXITER; niter++)
    {
        cout << niter << endl;
        ublas::vector<double> States_new = solver::TimeMarching(curved_mesh, param, resdata, States, invM, p, converged);
        States = States_new;
        if (converged)
            break;
    }
    return 0;
}
