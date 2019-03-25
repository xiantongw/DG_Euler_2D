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
#include "../include/CalcResidual.h"
#include "../include/euler.h"
#include "../include/Param.h"

using namespace std;
using namespace utils;
using namespace lagrange;


int main(int argc, char *argv[])
{
    namespace ublas = boost::numeric::ublas;

    string fname = "./mesh/bump0.gri";
    TriMesh mesh(fname);
    cout << mesh.gri_filename << endl;

    // cout << "Testing read file" << endl;
    // cout << mesh.E[mesh.num_element - 1][3] << ' ' << mesh.E[mesh.num_element - 1][4] << ' ' << mesh.E[mesh.num_element - 1][5] << endl;
    // cout << mesh.num_element << endl;
    // cout << mesh.name_base[0] << endl;
    // cout << endl;

    // cout << "Testing util" << endl;
    // int tmp[] = {1, 2, 3, 4, 5, 6};
    // std::vector<int> element(tmp, tmp + 6);
    // std::vector<int> test;
    // test = utils::GetVertexIndex(element);
    // cout << test[0] << ' ' << test[1] << ' ' << test[2] << endl;

    // int cl = 0;

    // cout << "Testing I2E" << endl;
    // cout << mesh.I2E[cl][0] << ' ' << mesh.I2E[cl][1] << ' ' << mesh.I2E[cl][2] << ' ' << mesh.I2E[cl][3] << ' ' << endl;
    // cout << endl;

    // cout << "Testing B2E" << endl;
    // cout << mesh.B2E[cl][0] << ' ' << mesh.B2E[cl][1] << ' ' << mesh.B2E[cl][2] << endl;
    // cout << endl;

    // for (int i = 0; i < 10; i++)
    // {
    //     for (int j = 0; j < 3; j++)
    //     {
    //         cout << mesh.B2E[i][j] << ' ';
    //     }
    //     cout << endl;
    // }

    // cout << endl
    //      << "Testing Lagrange" << endl;

    ublas::matrix<double> C;
    C = TriangleLagrange2D(2);
    // cout << C << endl;

    ublas::vector<ublas::vector<double> > vertex(3, ublas::vector<double>(2, 1));
    ublas::vector<ublas::vector<double> > physical(3, ublas::vector<double>(2, 1));
    ublas::vector<double> reference(2, 1), point(2, 1);

    vertex[0][0] = 0.0; vertex[0][1] = 0.0;
    vertex[1][0] = 1.0; vertex[1][1] = 0.0;
    vertex[2][0] = 0.0; vertex[2][1] = 1.0;

    point[0] = 0.5; point[1] = 0.5;

    physical = lagrange::MapReferenceToPhysicalLinear(vertex, 2);
    reference = lagrange::MapPhysicalToReferenceLinear(vertex, point, 2);

    // cout << vertex << endl;
    // cout << physical << endl;
    // cout << reference << endl;




    cout << endl <<  "You have entered " << argc - 1
         << " arguments:"
         << "\n";

    for (int i = 1; i < argc; ++i)
        cout << argv[i] << "\n";

    ublas::matrix<double> mat_mass = lagrange::ConstructMassMatrix(2, mesh);

    string boundary_name="bottom";

    // Testing Calculate Residaul
    cout << "Testing Calculate Residual" << endl;
    int p = 1;
    int Np = int((p + 1) * (p + 2) / 2);
    TriMesh curved_mesh = mesh;
    ConstructCurveMesh(mesh, curved_mesh, geometry::BumpFunction, boundary_name, p + 1);
    ublas::vector<double> States (curved_mesh.num_element * Np * 4, 0.0);
    // Construct free stream state
    for (int ielem = 0; ielem < curved_mesh.num_element; ielem++)
    {
        for (int ip = 0; ip < Np; ip++)
        {
            States(ielem * Np * 4 + ip * 4 + 0) = 1.0;
            States(ielem * Np * 4 + ip * 4 + 1) = 0.5;
            States(ielem * Np * 4 + ip * 4 + 2) = 0.0;
            States(ielem * Np * 4 + ip * 4 + 3) = 1.910714;
        }
    }
    Param param;
    // Set up the param struct
    /*************************/
    param.gamma = 1.40;
    param.attack_angle = 0.0;
    param.cfl = 0.5;
    param.mach_inf = 0.5;
    param.p_inf = 1.0;
    param.bound0 = "Free_Stream";
    param.bound1 = "Free_Stream";
    param.bound2 = "Free_Stream";
    param.bound3 = "Free_Stream";
    /*************************/
    ResData resdata = CalcResData(curved_mesh, p);
    ublas::vector<double> Residual = CalcResidual(curved_mesh, param, resdata, States, p);

    int ielem = 71;
    ublas::vector<double> elem_check(4, 0.0);
    int i_lagrange = 0;
    for (int istate = 0; istate < 4; istate++)
    {
        elem_check(istate) = Residual(ielem * Np * 4 + i_lagrange  * 4 + istate);
    }
    cout << "Total" << ublas::norm_2(elem_check) << endl;
    for (int i=0;i<curved_mesh.CurvedIndex.size();i++)
    {
        cout << curved_mesh.CurvedIndex[i] << ' ';
    }
    cout << endl;
    return 0;
}
