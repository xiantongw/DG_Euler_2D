#include "../include/lagrange.h"

namespace lagrange
{

    namespace ublas = boost::numeric::ublas;

    ublas::matrix<double> TriangleLagrange2D(int p)
    {
        /*
            calculates coeffs for full-order Lagrange basis of order preference element is a unit isoceles right triangle
        */
        ublas::vector<double> xi(p + 1);
        ublas::vector<double> eta(p + 1);
        double d = 1.0 / p;
        for (int i = 0; i < xi.size(); i++)
        {
            xi(i) = i * d;
            eta(i) = i * d;
        }
        int N = (p + 1) * (p + 2) / 2; // number of basis functions
        ublas::matrix<double> A(N, N);
        ublas::matrix<double> C(N, N);
        int i = 0, k = 0; // build A-matrix
        for (int iy = 0; iy <= p; iy++)
        {
            for (int ix = 0; ix <= p - iy; ix++)
            {
                k = 0;
                for (int s = 0; s <= p; s++)
                {
                    for (int r = 0; r <= p - s; r++)
                    {
                        A(i, k) = std::pow(xi(ix), r) * std::pow(eta(iy), s);
                        k = k + 1;
                    }
                }
                i = i + 1;
            }
        }
        bool inverted;
        inverted = InvertMatrix(A, C);
        if (!inverted)
        {
            std::cout << "Matrix cannot be inverted! Aborting..." << std::endl;
            abort();
        }
        return C;
    }


    ublas::vector<ublas::vector<double> > MapReferenceToPhysicalLinear(ublas::vector<ublas::vector<double> > vertex, int p)
    {
        int Np = int((p + 1) * (p + 2) / 2); // number of basis functions
        int r, s, k;
        ublas::vector<ublas::vector<double> > node_physical(Np, ublas::vector<double> (2, 1)), node_reference(Np, ublas::vector<double> (2, 1));
        ublas::matrix<double> jacobian(2, 2);

        // Initialize the Jacobian Matrix
        jacobian(0, 0) = vertex[1][0] - vertex[0][0];
        jacobian(0, 1) = vertex[2][0] - vertex[0][0];
        jacobian(1, 0) = vertex[1][1] - vertex[0][1];
        jacobian(1, 1) = vertex[2][1] - vertex[0][1];

        for (r = 0; r <= p; r++)
        {
            for (s = 0; s < p - r + 1; s++)
            {
                // Initialize the nodes in reference space
                k = utils::GetFullOrderIndex(r, s, p);
                node_reference[k][0] = 1.0 * r / p;
                node_reference[k][1] = 1.0 * s / p;
                // Map the nodes to the physical space
                /*  ->    ->   -   ->
                    x  = x_1 + J * xi*/
                ublas::axpy_prod(jacobian, node_reference[k], node_physical[k], true);
                node_physical[k] += vertex[0];
            }
        }
        return node_physical;
    }

    ublas::vector<double> MapPhysicalToReferenceLinear(ublas::vector<ublas::vector<double> > vertex,
                                                        ublas::vector<double> point, int p)
    {
        ublas::vector<double> node_reference(2, 1);
        ublas::matrix<double> jacobian(2, 2);
        ublas::matrix<double> inverse_jacobian(2, 2);
        double det_jacobian;

        // Initialize the Jacobian Matrix
        jacobian(0, 0) = vertex[1][0] - vertex[0][0];
        jacobian(0, 1) = vertex[2][0] - vertex[0][0];
        jacobian(1, 0) = vertex[1][1] - vertex[0][1];
        jacobian(1, 1) = vertex[2][1] - vertex[0][1];

        det_jacobian = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);

        inverse_jacobian(0, 0) = (vertex[2][1] - vertex[0][1]) / det_jacobian;
        inverse_jacobian(0, 1) = (vertex[0][0] - vertex[2][0]) / det_jacobian;
        inverse_jacobian(1, 0) = (vertex[0][1] - vertex[1][1]) / det_jacobian;
        inverse_jacobian(1, 1) = (vertex[1][0] - vertex[0][0]) / det_jacobian;

        ublas::axpy_prod(inverse_jacobian, point, node_reference, true);

        return node_reference;
    }

} // namespace lagrange
