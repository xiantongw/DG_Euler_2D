#include "../include/geometry.h"

namespace geometry
{
    double BumpFunction(double x)
    {
        return 0.0625 * std::exp(-25.0 * x * x);
    }

    ublas::matrix<double> CalcJacobianLinear(TriMesh mesh, int ielem)
    {
        ublas::matrix<double> jacobian (2, 2);
        ublas::matrix<double> vertex (3, 2);
        std::vector<int> vertex_index = utils::GetVertexIndex(mesh.E[ielem]);
        for (int i = 0; i < 3; i++)
        {
            vertex(i, 0) = mesh.V[mesh.E[ielem][vertex_index[i]] - 1][0];
            vertex(i, 1) = mesh.V[mesh.E[ielem][vertex_index[i]] - 1][1];
        }
        jacobian(0, 0) = vertex(1, 0) - vertex(0, 0);
        jacobian(0, 1) = vertex(2, 0) - vertex(0, 0);
        jacobian(1, 0) = vertex(1, 1) - vertex(0, 1);
        jacobian(1, 1) = vertex(2, 1) - vertex(0, 1);

        return jacobian;
    }

    ublas::matrix<double> CalcJacobianCurved(TriMesh mesh, int ielem, boost::multi_array<double, 3> GPhi_Curved, int n_quad_2d, int ig)
    {
        // In a curved element, the jacobian varies from the points eveluated on
        // Since GPhi on quadrature points are pre-caculated, so GPhi is passed in
        ublas::matrix<double> jacobian(2, 2, 0.0);
        // First, if this element is really a curved element?
        std:vector<int> lagrange_nodes_index = mesh.E[ielem];
        int Np = lagrange_nodes_index.size();
        int p = int((sqrt(1 + 8.0 * Np) - 3) / 2);
        ublas::matrix<double> Nodes_Coord(Np, 2);
        for (int i = 0; i < Np; i++)
        {
            Nodes_Coord(i, 0) = mesh.V[lagrange_nodes_index[i] - 1][0];
            Nodes_Coord(i, 1) = mesh.V[lagrange_nodes_index[i] - 1][1];
        }
        ublas::matrix <double> GPhi_on_quad(Np, 2, 0.0);
        GPhi_Curved.resize(boost::extents[n_quad_2d][Np][2]);
        for (int ip = 0; ip < Np; ip++)
        {
            GPhi_on_quad(ip, 0) = GPhi_Curved[ig][ip][0];
            GPhi_on_quad(ip, 1) = GPhi_Curved[ig][ip][1];
        }
        ublas::axpy_prod(ublas::trans(Nodes_Coord), GPhi_on_quad, jacobian, true);
        return jacobian;
    }

    ublas::vector<int> GetEdgeLagrangeNodeIndex(int p, int local_index)
    {
        ublas::vector<int> selected_lagrange_index(p - 1);
        switch (local_index)
        {
            case 0:
                if (p > 1)
                {
                    selected_lagrange_index[0] = p + 1 + p + 1 - 1;
                    for (int i = 1; i < p - 1; i++)
                    {
                        selected_lagrange_index[i] = selected_lagrange_index[i - 1] + p + 1 - i - 1;
                    }
                }
                for (int i = 0; i < selected_lagrange_index.size(); i++)
                {
                    selected_lagrange_index[i] = selected_lagrange_index[i] - 1;
                }
                break;

            case 1:
                if (p > 1)
                {
                    selected_lagrange_index[0] = p + 1;
                    for (int i = 1; i < p - 1; i++)
                    {
                        selected_lagrange_index[i] = selected_lagrange_index[i - 1] + 1 + p - i;
                    }
                }
                boost::range::reverse(selected_lagrange_index);
                break;

            case 2:
                for (int i = 0; i < p - 1; i++)
                {
                    selected_lagrange_index[i] = i + 1;
                }
                break;

            default:
                cout << "Mesh Error! Aborting..." << endl;
                abort();
        }
        return selected_lagrange_index;
    }

    ublas::vector<int> GetEdgeCoordinatesIndex(TriMesh mesh, int iedge)
    {
        int Nq = mesh.E[mesh.B2E[iedge][0] - 1].size();
        int q = int((sqrt(1 + 8.0 * Nq) - 3) / 2);
        int ielemL = mesh.B2E[iedge][0] - 1;
        int ilocL = mesh.B2E[iedge][1] - 1;
        ublas::vector<int> selected_lagrange_index = GetEdgeLagrangeNodeIndex(q, ilocL);
        ublas::vector<int> edge_coord_index(selected_lagrange_index.size() + 2, 0);
        int iL, iR;
        switch (ilocL)
        {
            case 0:
                iL = q; iR = Nq - 1;
                break;
            case 1:
                iL = Nq - 1; iR = 0;
                break;
            case 2:
                iL = 0; iR = q;
            default:
                break;
        }
        edge_coord_index(0) = iL;
        edge_coord_index(selected_lagrange_index.size() + 1) = iR;
        for (int i = 1; i < selected_lagrange_index.size() + 1; i++)
        {
            edge_coord_index(i) = selected_lagrange_index[i - 1];
        }
        return edge_coord_index;
    }

    ublas::vector<ublas::vector<double> > GetEdgeCoordinates(TriMesh mesh, int iedge)
    {
        // Get the geometry points on the edge
        int Nq = mesh.E[mesh.B2E[iedge][0] - 1].size();
        int q = int((sqrt(1 + 8.0 * Nq) - 3) / 2);
        int ielemL = mesh.B2E[iedge][0] - 1;
        int ilocL = mesh.B2E[iedge][1] - 1;
        ublas::vector<ublas::vector<double> > edge_coord(q + 1, ublas::vector<double> (2, 0.0));
        ublas::vector<int> edge_coord_index = GetEdgeCoordinatesIndex(mesh, iedge);

        for (int i = 0; i < edge_coord_index.size(); i++)
        {
            edge_coord(i)(0) = mesh.V[mesh.E[ielemL][edge_coord_index[i]] - 1][0];
            edge_coord(i)(1) = mesh.V[mesh.E[ielemL][edge_coord_index[i]] - 1][1];
        }
        return edge_coord;
    }

} // namespace geometry
