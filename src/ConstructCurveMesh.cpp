#include "../include/ConstructCurveMesh.h"

TriMesh ConstructCurveMesh(TriMesh& mesh, TriMesh& curved_mesh, double (*pBumpFunction)(double), string boundary_name, int p)
{
    namespace ublas = boost::numeric::ublas;
    int boundary_index;
    // Find the index of the boundary which need to be curved
    int match = 0;
    for(int i = 0; i < mesh.Bname.size(); i++)
    {
        if (strcasecmp(mesh.Bname[i].c_str(), boundary_name.c_str()) == 0)
        {
            boundary_index = i;
            match = 1;
            break;
        }
    }
    if (!match)
    {
        cout << "Undefined boundary name! " << mesh.Bname[0].c_str() << " Aborting..." << endl;
        abort();
    }

    int num_boundary_elements = mesh.B[boundary_index].size();
    // construct the ublas vector based on the std vector in the TriMesh
    ublas::vector<ublas::vector<int> > boundary_elements(num_boundary_elements, ublas::vector<int> (3));
    int element_counter_boundary = 0, element_counter_total = 0;
    while (element_counter_boundary < num_boundary_elements)
    {
        if (mesh.B2E[element_counter_total][2] == boundary_index + 1)
        {
            for (int i = 0; i < 3; i++)
            {
                boundary_elements[element_counter_boundary][i] = mesh.B2E[element_counter_total][i];
            }
            element_counter_boundary++;
        }
        element_counter_total++;
    }
    std::vector<int> index_boundary_elements (num_boundary_elements);
    for (int i = 0; i < num_boundary_elements; i++)
    {
        index_boundary_elements[i] = boundary_elements[i][0] - 1;
    }
    // The index of three vetex, which shouldn't be changed
    ublas::vector<int> index_kept(3);
    index_kept[0] = utils::GetFullOrderIndex(0, 0, p);
    index_kept[1] = utils::GetFullOrderIndex(p, 0, p);
    index_kept[2] = utils::GetFullOrderIndex(0, p, p);

    // Loop through all elements and do curve on boundary elements
    for (int i_elem = 0; i_elem < mesh.E.size(); i_elem++)
    {
        ublas::vector<int> v_ind (3); // indices of the vertex
        ublas::vector<ublas::vector<double> > node_lagrange, vertex(3, ublas::vector<double> (2));

        for (int i = 0; i < 3; i++)
        {
            v_ind[i] = mesh.E[i_elem][i] - 1;
            vertex[i][0] = mesh.V[v_ind[i]][0];
            vertex[i][1] = mesh.V[v_ind[i]][1];
        }
        node_lagrange = lagrange::MapReferenceToPhysicalLinear(vertex, p);

        // Check if the element need to be curved
        std::vector<int>::iterator it = std::find(index_boundary_elements.begin(), index_boundary_elements.end(), i_elem);
        if(it != index_boundary_elements.end())
        {
            // First, find which nodes in node_lagrange are on the boundary need to be curved
            int index_in_boudary_elements = std::distance(index_boundary_elements.begin(), it);
            int local_index = boundary_elements[index_in_boudary_elements][1] - 1;
            ublas::vector<int> selected_lagrange_index = geometry::GetEdgeLagrangeNodeIndex(p, local_index);
            // modify the y coordinate using the bump function
            for (int i = 0; i < selected_lagrange_index.size(); i++)
            {
                node_lagrange[selected_lagrange_index[i]][1] = pBumpFunction(node_lagrange[selected_lagrange_index[i]][0]);
            }
            curved_mesh.isCurved[i_elem] = true;

            // Now we have node_lagrange, which is the new nodes need to be written into the curved_mesh.B and we add new points
            // NOTICE: The index of 3 vertex shouldn't be changed
            curved_mesh.E[i_elem] = std::vector<int> (int((p + 1) * (p + 2) / 2), 0);

            for (int i_kept = 0; i_kept < 3; i_kept++)
            {
                curved_mesh.E[i_elem][index_kept[i_kept]] = mesh.E[i_elem][i_kept];
            }
            // add points to curved_mesh.V, !!! NOT the current vertex !!!
            for (int i = 0; i < node_lagrange.size(); i++)
            {
                if ((i != index_kept[0]) && (i != index_kept[1]) && (i != index_kept[2]))
                {
                    std::vector<double> new_point {node_lagrange[i][0], node_lagrange[i][1]};
                    // find if the new point alreay exists
                    int j = 0;
                    for (j = 0; j < curved_mesh.V.size(); j++)
                    {
                        if (std::abs(new_point[0] - curved_mesh.V[j][0]) < 1e-6 && std::abs(new_point[1] - curved_mesh.V[j][1]) < 1e-6)
                        {
                            curved_mesh.E[i_elem][i] = j + 1;
                            break;
                        }
                    }
                    if (j == curved_mesh.V.size())
                    {
                        curved_mesh.V.push_back(new_point);
                        curved_mesh.E[i_elem][i] = curved_mesh.V.size();
                    }
                }
            }
        } // end of dealing with curved element
    } // End loop through all elements
    curved_mesh.curved_group = boundary_index + 1;
    curved_mesh.FindCurvedIndex(); // Calculate the index of the curved elements
    curved_mesh.CalcBn();
    return curved_mesh;
}
