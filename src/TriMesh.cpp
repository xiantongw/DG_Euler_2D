#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "../include/utils.h"
#include "../include/TriMesh.h"

TriMesh::TriMesh(string& gri_filename_in){
    gri_filename = gri_filename_in;
    ReadGri(this->gri_filename);
    calculateI2E();
    calculateB2E();
}

void TriMesh::ReadGri(string& gri_filename){

    // variables dealing with strings
    ifstream gri_file(gri_filename);
    stringstream ss;
    string line;
    // Variables
    int num_node, num_element, dim, num_boundary, n_base;
    // Read the first line, get Nn, Ne, dim
    getline(gri_file, line); ss.clear(); ss.str(line);
    ss >> num_node >> num_element >> dim;

    // Read the V matrix, which are the coordinates of the nodes
    vector<vector<double> > V(num_node, vector<double>(2));
    for (int i = 0; i < num_node; i++){
        getline(gri_file, line); ss.clear(); ss.str(line);
        ss >> V[i][0] >> V[i][1];
    }

    // Read number of boundaries
    getline(gri_file, line); ss.clear(); ss.str(line);
    ss >> num_boundary;

    // Read boundary information
    vector<string> Bname(num_boundary);
    vector<vector<vector<int> > > B;
    B.resize(num_boundary);
    for (int i_boundary = 0; i_boundary < num_boundary; i_boundary++){
        // The information of a specific boundary 
        int num_boundary_edge, dim;
        getline(gri_file, line); ss.clear(); ss.str(line);
        ss >> num_boundary_edge >> dim >> Bname[i_boundary];
        B[i_boundary].resize(num_boundary_edge);
        for (int i = 0; i < num_boundary_edge; i++){
            getline(gri_file, line); ss.clear(); ss.str(line);
            B[i_boundary][i].resize(2);
            ss >> B[i_boundary][i][0] >> B[i_boundary][i][1];
        }
    }

    // Read element information
    vector<vector<int> > E(num_element);
    vector<string> name_base;
    vector<int> n_bases;

    int current_total_elements = 0, num_element_part, num_element_part_old = 0;
    int nnode, buffer_int;
    string buffer_str;

    while (current_total_elements != num_element){

        getline(gri_file, line); ss.clear(); ss.str(line);
        ss >> num_element_part >> n_base >> buffer_str;
        n_bases.push_back(n_base);
        name_base.push_back(buffer_str);

        for (int i_element = 0; i_element < num_element_part; i_element++){
            getline(gri_file, line); ss.clear(); ss.str(line);
            nnode = (n_base + 1) * (n_base + 2) / 2;
            
            E[i_element + num_element_part_old].resize(nnode);
            int i_local_node = 0; 
            while (ss >> buffer_int){
                E[i_element + num_element_part_old][i_local_node] = buffer_int;
                i_local_node ++;
            }
        }
        current_total_elements += num_element_part;
        num_element_part_old = num_element_part;
    }
    gri_file.close();

    this->B = B;
    this->Bname = Bname;
    this->E = E;
    this->V = V;
    this->num_boundary = num_boundary;
    this->num_element = num_element;
    this->num_node = num_node;
    this->name_base = name_base;
    this->n_base = n_bases;

}

void TriMesh::calculateI2E(){
    // This function is re-written from the Python version, so the names of the 
    // variables are not well-organized
    vector<vector<int> > E = this->E;
    vector<vector<double> > V = this->V;
    int N = this->num_node;
    boost::numeric::ublas::mapped_matrix<int> H(N, N, N * N);
    vector<vector<int> > C(this->num_element * 3, vector<int>(4, 0));
    int nedge = 0;
    for (int t = 0; t < this->num_element; t++){
        vector<int> vertex_index(3);
        vertex_index = utils::getVertexIndex(E[t]);
        for (int e = 0; e < 3; e++){
            int n1 = E[t][vertex_index[(e + 1) % 3]] - 1; 
            int n2 = E[t][vertex_index[(e + 2) % 3]] - 1;
            int nmin = min(n1, n2); int nmax = max(n1, n2);
            if (H(nmin, nmax) > 0){
                int tN = H(nmin, nmax); int eN = H(nmax, nmin);
                int t1 = t + 1; int t2 = tN;
                int e1 = e + 1;  int e2 = eN;
                if (t2 < t1) {
                    t1 = tN; t2 = t + 1;
                    e1 = eN; e2 = e + 1; 
                }
                C[nedge][0] = t1; C[nedge][1] = e1;
                C[nedge][2] = t2; C[nedge][3] = e2;
                nedge = nedge + 1;
            }
            else {
                H(nmin, nmax) = t + 1;
                H(nmax, nmin) = e + 1;
            }
        }
    }
    // Sorting
    vector<vector<int> > CC(nedge, vector<int>(4, 0));
    for (int i = 0; i < nedge; i++){
        for (int j = 0; j < 4; j++){
            CC[i][j] = C[i][j];
        }
    }
    sort(CC.begin(), CC.end(), utils::sortcol_0);
    int i0 = 0;
    for (int i = 0; i < nedge; i++){
        if (CC[i][0] != CC[i0][0]){
            vector<vector<int> > A = utils::slice_by_row(CC, i0, i - 1);
            sort(A.begin(), A.end(), utils::sortcol_2);
            for (int ii = i0; ii < i; ii++){
                for (int jj = 0; jj < 4; jj++){
                    CC[ii][jj] = A[ii - i0][jj];
                }
            }
            i0 = i;
        }
    }
    this->I2E = CC;
}


void TriMesh::calculateB2E(){
    // This function is re-written from the Python version, so the names of the 
    // variables are not well-organized
    vector<vector<int> > E = this->E;
    vector<vector<vector<int> > > B = this->B;
    vector<vector<int> > B2E;
    
    for (int ielem = 0; ielem < E.size(); ielem++){
        vector<int> vertex_index(3), elem(3);
        vertex_index = utils::getVertexIndex(E[ielem]);
        for (int i = 0; i < 3; i++){
            elem[i] = E[ielem][i];
        }
        for (int ibg = 0; ibg < this->Bname.size(); ibg++){
            vector<vector<int> > B_part = B[ibg];
            for (int ib = 0; ib < B_part.size(); ib++){
                vector<int> nb = B_part[ib];
                int match_1 = distance(elem.begin(), find(elem.begin(), elem.end(), nb[0]));
                int match_2 = distance(elem.begin(), find(elem.begin(), elem.end(), nb[1]));
                if ((match_1 != elem.size()) && (match_2 != elem.size())){
                    int local_ind = utils::missing_from_012(match_1, match_2);
                    vector<int> row = {ielem + 1, local_ind + 1, ibg + 1};
                    B2E.push_back(row);
                }
            }
        }
    }
    sort(B2E.begin(), B2E.end(), utils::sortcol_2);
    this->B2E = B2E;
}
