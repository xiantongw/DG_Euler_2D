#ifndef TRIMESH_H
#define TRIMESH_H

#include <string>
#include <vector>

using namespace std;

class TriMesh
{
    public:
    
        string gri_filename;
        double num_node;
        int num_element;
        double num_boundary;
        vector<vector<int> > E;
        vector<vector<double> > V;
        vector<vector<vector<int> > > B;
        vector<string> Bname;
        vector<string> name_base;
        vector<int> n_base;
        vector<vector<int> > I2E;
        vector<vector<int> > B2E;
        vector<vector<double> > In;
        vector<vector<double> > Bn;
        
        TriMesh(string& gri_filename_in);
        void ReadGri(string& gri_filename);
        void CalculateI2E();
        void CalculateB2E();

};

#endif