#include <iostream>
#include <fstream>
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

int main(int argc, char *argv[])
{
    Param param;
    // Set up the param struct
    param = ReadParamIn(string(argv[1]));
    TriMesh mesh(param.mesh_file);
    // Testing Calculate Residaul
    int p = param.order;
    int q = param.order_geo;
    int Np = int((p + 1) * (p + 2) / 2);
    TriMesh curved_mesh = mesh;
    string boundary_name="bottom";
    ConstructCurveMesh(mesh, curved_mesh, geometry::BumpFunction, boundary_name, q);

    ResData resdata_postproc;
    if (p == 0)
    {
        solver::CalcResData(curved_mesh, 1, resdata_postproc);
    }else
    {
        solver::CalcResData(curved_mesh, p, resdata_postproc);
    }

    int Np_solution;
    if (p == 0)
        Np_solution = 3;
    else
        Np_solution = Np;
    ublas::vector<ublas::matrix<double> > Nodes(curved_mesh.E.size(), ublas::matrix<double>(Np_solution, 2, 0.0));
    ublas::vector<ublas::matrix<double> > State_on_Nodes(curved_mesh.E.size(), ublas::matrix<double>(Np_solution, 4, 0.0));

    ifstream file_nodes, file_state;
    file_nodes.open("nodes.dat");
    file_state.open("states.dat");
    for (int ielem = 0; ielem < curved_mesh.E.size(); ielem++)
    {
        for (int ip = 0; ip < Np_solution; ip++)
        {
            double x, y, s0, s1, s2, s3;
            char temp;
            file_nodes  >> x  >> y;
            file_state  >> s0 >> s1 >> s2 >> s3;
            Nodes(ielem)(ip, 0) = x; Nodes(ielem)(ip, 1) = y;
            State_on_Nodes(ielem)(ip, 0) = s0; State_on_Nodes(ielem)(ip, 1) = s1;
            State_on_Nodes(ielem)(ip, 2) = s2; State_on_Nodes(ielem)(ip, 3) = s3;
        }
    }
    file_nodes.close();
    file_state.close();
    param.h = 0.0625;
    double err_entropy, coeff_lift, coeff_drag;
    std::vector<std::vector<double> > p_coeff_dist;
    solver::CalcScalarOutputs(curved_mesh, resdata_postproc, State_on_Nodes, Nodes, param,
                                    err_entropy, coeff_lift, coeff_drag, p_coeff_dist);
    cout.setf(ios::scientific, ios::floatfield);
    std::cout << setprecision(8) << err_entropy << ' ' << coeff_lift << ' ' << coeff_drag << endl;
    ofstream file_pdist;
    file_pdist.open("pressure_distribution.dat");
    for (int i = 0; i < p_coeff_dist.size(); i++)
    {
        cout.setf(ios::scientific, ios::floatfield);
        file_pdist  << setprecision(10) << p_coeff_dist[i][0] << ' ' << p_coeff_dist[i][1] << std::endl;
    }
    file_pdist.close();
    return 0;
}
