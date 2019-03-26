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
        return ublas::trans(C);
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

    ublas::vector<ublas::matrix<double> > ConstructMassMatrix(int p, TriMesh mesh, ResData resdata)
    {
        int Np = int((p + 1) * (p + 2) / 2);
        int num_element = mesh.E.size();
        ublas::vector<ublas::matrix<double> > mat_mass (num_element, ublas::matrix<double> (Np, Np, 0.0));
        ublas::matrix<double> unit_mat_mass (Np, Np, 0.0);
        std::vector<double> xq_std;
        std::vector<double> wq_std;
        int n_quad_2d;
        GetQuadraturePointsWeight2D(2, n_quad_2d, xq_std, wq_std);
        // transfer the std vectors to ublas vectors
        ublas::vector<double> xq = utils::StdToBoostVector(xq_std);
        ublas::vector<double> wq = utils::StdToBoostVector(wq_std);
        ublas::matrix<double> TriLagrangeCoeff = TriangleLagrange2D(p);
        // Compute the base functions on all quadrature points
        ublas::matrix<double> Phi(n_quad_2d, Np); // Basis function value on qudatrue points
        for (int ig = 0; ig < n_quad_2d; ig++)
        {
            double xi = xq(2 * ig);
            double eta = xq(2 * ig + 1);
            ublas::vector<double> phi = CalcBaseFunction(TriLagrangeCoeff, xi, eta);
            for (int ip = 0; ip < Np; ip++)
            {
                Phi(ig, ip) = phi(ip);
            }
        }
        // Do the Gussian quadrature integration for linear element
        unit_mat_mass.clear();
        for (int i = 0; i < Np; i++)
        {
            for (int j = 0; j < Np; j++)
            {
                for (int ig = 0; ig < n_quad_2d; ig++)
                {
                    // unit mass matrix is the same in reference space is the same for linear elements
                    unit_mat_mass(i, j) += Phi(ig, i) * Phi(ig, j);
                }
            }
        }
        // Fill the term of block-diagonal parts for linear elements
        for (int i_elem = 0; i_elem < num_element; i_elem++)
        {
            if (!mesh.isCurved[i_elem])
            {
                for (int i = 0; i < Np; i++)
                {
                    for (int j = 0; j < Np; j++)
                    {
                        ublas::matrix<double> jacobian = geometry::CalcJacobianLinear(mesh, i_elem);
                        double det_jacobian = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
                        mat_mass (i_elem)(i, j) = unit_mat_mass (i, j) * det_jacobian;
                    }
                }
            }
        }
        // Fill the term of block-diagonal parts for curved elements
        for (int i_elem = 0; i_elem < num_element; i_elem++)
        {
            if (mesh.isCurved[i_elem])
            {
                for (int i = 0; i < Np; i++)
                {
                    for (int j = 0; j < Np; j++)
                    {
                        unit_mat_mass.clear();
                        double det_jacobian = 0.0;
                        for (int ig = 0; ig < n_quad_2d; ig++)
                        {
                            // unit mass matrix is the same in reference space is the same for linear elements
                            ublas::matrix<double> jacobian = geometry::CalcJacobianCurved(mesh, i_elem, resdata.GPhi_Curved, n_quad_2d, ig);
                            det_jacobian = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
                            unit_mat_mass(i, j) += Phi(ig, i) * Phi(ig, j) * det_jacobian;
                        }
                        mat_mass (i_elem)(i, j) = unit_mat_mass (i, j);
                    }
                }
            }
        }
        return mat_mass;
    }

    ublas::vector<ublas::matrix<double> > CalcInvMassMatrix(ublas::vector<ublas::matrix<double> > M)
    {
        int num_elem = M.size();
        int Np = M(0).size1();
        ublas::vector<ublas::matrix<double> > invM(num_elem, ublas::matrix<double>(Np, Np, 0.0));
        for (int ielem = 0; ielem < num_elem; ielem++)
        {
            ublas::matrix<double> unit_invM (Np, Np, 0.0);
            ublas::matrix<double> unit_M = M(ielem);
            InvertMatrix(unit_M, unit_invM);
            invM(ielem) = unit_invM;
        }
        return invM;
    }

    ublas::vector<double> CalcBaseFunction(ublas::matrix<double> TriLagrangeCoeff, double xi, double eta)
    {
        int num_poly = TriLagrangeCoeff.size2();
        int p = int((sqrt(1 + 8.0 * num_poly) - 3) / 2);
        int s, r, ind = 0;
        ublas::vector<double> monomial(num_poly);
        for (s = 0; s <= p; s++)
        {
            for (r = 0; r <= (p - s); r++)
            {
                monomial(ind) = pow(xi, r) * pow(eta, s);
                ind++;
            }
        }
        ublas::vector<double> phi(TriLagrangeCoeff.size1());
        for (int i = 0; i < TriLagrangeCoeff.size1(); i++)
        {
            double temp_sum = 0.0;
            for (int j = 0; j < TriLagrangeCoeff.size2(); j++)
            {
                temp_sum += TriLagrangeCoeff(i, j) * monomial(j);
            }
            phi(i) = temp_sum;
        }
        return phi;
    }

    ublas::matrix<double> CalcBaseFunctionGradient(ublas::matrix<double> TriLagrangeCoeff, double xi, double eta)
    {
        int num_poly = TriLagrangeCoeff.size2();
        int p = int((sqrt(1 + 8.0 * num_poly) - 3) / 2);
        int s, r, ind = 0;
        ublas::matrix<double> monomial(num_poly, 2, 0.0);
        for (s = 0; s <= p; s++)
        {
            for (r = 0; r <= (p - s); r++)
            {
                if (abs(xi) < 1e-6 && r == 0 && s > 0)
                {
                    monomial(ind, 0) = 0.0;
                    monomial(ind, 1) = s * pow(xi, r) * pow(eta, s - 1);
                }
                else if (abs(eta) < 1e-6 && s == 0 && r > 0)
                {
                    monomial(ind, 0) = r * pow(xi, r - 1) * pow(eta, s);
                    monomial(ind, 1) = 0.0;
                }
                else if (r == 0 && s == 0)
                {
                    monomial(ind, 0) = 0.0;
                    monomial(ind, 1) = 0.0;
                }
                else
                {
                    monomial(ind, 0) = r * pow(xi, r - 1) * pow(eta, s);
                    monomial(ind, 1) = s * pow(xi, r) * pow(eta, s - 1);
                }
                ind++;
            }
        }
        ublas::matrix<double> phi(TriLagrangeCoeff.size1(), 2);
        for (int i = 0; i < TriLagrangeCoeff.size1(); i++)
        {
            double temp_sum[2] = {0.0};
            for (int j = 0; j < TriLagrangeCoeff.size2(); j++)
            {
                temp_sum[0] += TriLagrangeCoeff(i, j) * monomial(j, 0);
                temp_sum[1] += TriLagrangeCoeff(i, j) * monomial(j, 1);
            }
            phi(i, 0) = temp_sum[0];
            phi(i, 1) = temp_sum[1];
        }
        return phi;
    }

} // namespace lagrange
