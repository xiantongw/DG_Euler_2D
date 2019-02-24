#include <iostream>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "../include/TriMesh.h"
#include "../include/utils.h"

using namespace std;
using namespace utils;

int main()
{

    string fname = "./mesh/bump2.gri";
    TriMesh mesh(fname);
    cout << mesh.gri_filename << endl;
    
    cout << "Testing read file" << endl;
    cout << mesh.E[mesh.num_element-1][0] << ' ' << mesh.E[mesh.num_element-1][1] << ' ' << mesh.E[mesh.num_element-1][2] <<endl;
    cout << mesh.num_element <<endl;
    cout << mesh.name_base[0] << endl;
    cout << endl;

    
    cout << "Testing util" << endl;
    int tmp[] = {1, 2, 3, 4, 5, 6};
    vector<int> element (tmp, tmp+6);
    vector<int> test;
    test = utils::GetVertexIndex(element);
    cout << test[0] << ' ' << test[1] << ' ' << test[2] << endl; 

    cout << endl << "Testing boost" << endl;
    boost::numeric::ublas::mapped_matrix<double> m (3, 3, 3 * 3);
    for (unsigned i = 0; i < m.size1 (); ++ i)
        for (unsigned j = 0; j < m.size2 (); ++ j)
            m (i, j) = 3 * i + j;
    cout << m << endl;

    
    cout << "Testing Slice" << endl;
    vector< vector<int> > vect(3, vector<int>(3)); 

	vect[0][0] = 3; vect[0][1] = 5; vect[0][2] = 1;
	vect[1][0] = 4; vect[1][1] = 8; vect[1][2] = 6;
	vect[2][0] = 7; vect[2][1] = 2; vect[2][2] = 9;

    cout << vect[0][0] << vect[0][1] << vect[0][2] << endl;
    
    
    vector<vector<int> > sub_vec = utils::SliceByRow(vect, 0, 1); 
    
    cout << sub_vec[0][0] << sub_vec[0][1] << sub_vec[0][2] << endl;
    cout << sub_vec[1][0] << sub_vec[1][1] << sub_vec[1][2] << endl;

    cout << endl;
    
    int cl = 0;
    
    cout << "Testing I2E" << endl;
    cout << mesh.I2E[cl][0] << ' ' << mesh.I2E[cl][1] << ' ' << mesh.I2E[cl][2] << ' ' << mesh.I2E[cl][3] << ' ' <<endl;
    cout << endl;

    cout << "Testing B2E" << endl;
    cout << mesh.B2E[cl][0] << ' ' << mesh.B2E[cl][1] << ' ' << mesh.B2E[cl][2] <<endl;
    cout << endl;

    for (int i = 0; i < 10; i++){
        for (int j=0; j<3; j++){
            cout << mesh.B2E[i][j] << ' ';
        }
        cout << endl;
    }
}
