#include <iostream>
#include <algorithm>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "../include/TriMesh.h"
#include "../include/utils.h"
#include "../include/lagrange.h"
#include "../include/geometry.h"
#include "../include/ConstructCurveMesh.h"

using namespace std;
using namespace utils;
using namespace lagrange;


int main(int argc, char *argv[])
{
    namespace ublas = boost::numeric::ublas;

    string fname = "./mesh/bump0.gri";
    TriMesh mesh(fname);
    cout << mesh.gri_filename << endl;

    cout << "Testing read file" << endl;
    cout << mesh.E[mesh.num_element - 1][0] << ' ' << mesh.E[mesh.num_element - 1][1] << ' ' << mesh.E[mesh.num_element - 1][2] << endl;
    cout << mesh.num_element << endl;
    cout << mesh.name_base[0] << endl;
    cout << endl;

    cout << "Testing util" << endl;
    int tmp[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> element(tmp, tmp + 6);
    std::vector<int> test;
    test = utils::GetVertexIndex(element);
    cout << test[0] << ' ' << test[1] << ' ' << test[2] << endl;

    cout << endl
         << "Testing boost" << endl;
    boost::numeric::ublas::mapped_matrix<double> m(3, 3, 3 * 3);
    for (unsigned i = 0; i < m.size1(); ++i)
        for (unsigned j = 0; j < m.size2(); ++j)
            m(i, j) = 3 * i + j;
    cout << m << endl;

    int cl = 0;

    cout << "Testing I2E" << endl;
    cout << mesh.I2E[cl][0] << ' ' << mesh.I2E[cl][1] << ' ' << mesh.I2E[cl][2] << ' ' << mesh.I2E[cl][3] << ' ' << endl;
    cout << endl;

    cout << "Testing B2E" << endl;
    cout << mesh.B2E[cl][0] << ' ' << mesh.B2E[cl][1] << ' ' << mesh.B2E[cl][2] << endl;
    cout << endl;

    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << mesh.B2E[i][j] << ' ';
        }
        cout << endl;
    }

    cout << endl
         << "Testing Lagrange" << endl;

    ublas::matrix<double> C;
    C = TriangleLagrange2D(2);
    cout << C << endl;

    ublas::vector<ublas::vector<double> > vertex(3, ublas::vector<double>(2, 1));
    ublas::vector<ublas::vector<double> > physical(3, ublas::vector<double>(2, 1));
    ublas::vector<double> reference(2, 1), point(2, 1);

    vertex[0][0] = 0.0; vertex[0][1] = 0.0;
    vertex[1][0] = 1.0; vertex[1][1] = 0.0;
    vertex[2][0] = 0.0; vertex[2][1] = 1.0;

    point[0] = 0.5; point[1] = 0.5;

    physical = lagrange::MapReferenceToPhysicalLinear(vertex, 2);
    reference = lagrange::MapPhysicalToReferenceLinear(vertex, point, 2);

    cout << vertex << endl;
    cout << physical << endl;
    cout << reference << endl;




    cout << endl <<  "You have entered " << argc - 1
         << " arguments:"
         << "\n";

    for (int i = 1; i < argc; ++i)
        cout << argv[i] << "\n";

    cout << endl << "Testing geometry" << endl;
    cout << geometry::BumpFunction(0.1) << endl;

    cout << endl << "test construct curve mesh" << endl;
    string boundary_name="bottom";
    cout << mesh.Bname[0] << endl;
    TriMesh curved_mesh = mesh;
    ConstructCurveMesh(mesh, curved_mesh, geometry::BumpFunction, boundary_name, 4);

    cout << mesh.V.size() << "||" << curved_mesh.V.size() << endl;
    cout << mesh.E.size() << "||" << curved_mesh.E.size() << endl;

    cout << endl;
    string out_filename = "test_out.gri";
    curved_mesh.WriteGri(out_filename);

    cout << mesh.isCurved[1] << endl;

    return 0;
}
