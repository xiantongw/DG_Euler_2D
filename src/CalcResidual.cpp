#include "../include/CalcResidual.h"

namespace ublas = boost::numeric::ublas;

ublas::vector<double> CalcResidual(TriMesh mesh, Param& param, ResData& resdata, ublas::vector<double> States, int p)
{
    // Unroll the mesh information
    int num_element = mesh.num_element;
    int num_states = 4; // Four states for 2D Euler equation
    // The number of lagrange nodes in each element
    int Np = int((p + 1) * (p + 2) / 2);
    int Nq = mesh.E[mesh.CurvedIndex[0]].size();
    int q = int((sqrt(1 + 8.0 * Nq) - 3) / 2);
    double gamma = param.gamma;
    // The residual vector
    ublas::vector<double> Residual (num_element * Np * num_states, 0.0);

    // Pull out the ResData
    arr_2d Phi = resdata.Phi;
    arr_3d GPhi = resdata.GPhi;
	arr_3d Phi_1D = resdata.Phi_1D;
	arr_3d Phi_1D_Curved = resdata.Phi_1D_Curved;
	arr_4d GPhi_1D = resdata.GPhi_1D;
	arr_4d GPhi_1D_Curved = resdata.GPhi_1D_Curved;
	int n_quad_1d = resdata.n_quad_1d;
	int n_quad_2d = resdata.n_quad_2d;
	ublas::vector<double> x_quad_1d = resdata.x_quad_1d;
	ublas::vector<double> w_quad_1d = resdata.w_quad_1d;
	ublas::vector<double> x_quad_2d = resdata.x_quad_2d;
	ublas::vector<double> w_quad_2d = resdata.w_quad_2d;

    // Loop over elements
    // 1. Construct Jacobian 2.evealute interior contribution to the residual
    // Notice that there are two sets of elements in the mesh, curved and not curved
    for (int ielem = 0; ielem < num_element; ielem++)
    {
        // pull out the state from the big State vector on each lagrange node which will be used on
        // interpolation to quadrature points
        ublas::vector<ublas::vector<double> > states_in_element(Np, ublas::vector<double> (num_states, 0.0));
        for (int ip = 0; ip < Np; ip++)
        {
            ublas::vector<double> state(num_states);
            for (int istate = 0; istate < num_states; istate++)
            {
                state(istate) = States(ielem * Np * num_states + ip * num_states + istate);
            }
            states_in_element(ip) = state;
        }
        if (! mesh.isCurved[ielem]) // Not a curved element
        {
            /*********************************************/
            /* Interior Contribution from linear element */
            /*********************************************/
            // Construct the Jacobian matrix
            ublas::matrix<double> jacobian = geometry::CalcJacobianLinear(mesh, ielem);
            ublas::matrix<double> inv_jacobian (jacobian.size1(), jacobian.size2());
            double det_jacobian = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
            InvertMatrix(jacobian, inv_jacobian);
            // Do the integration using quadrature points, which is the interior contribution of the residual
            for (int ip = 0; ip < Np; ip++) // loop over all lagrange nodes
            {
                ublas::vector<double> temp_sum (num_states, 0.0); // the temp summation for the quadrature integration
                for (int ig = 0; ig < n_quad_2d; ig++) // loop over all quadrature points
                {
                    // interpolate the state to quadrature points
                    ublas::vector<double> state_on_quad(num_states, 0.0);
                    for (int ipi = 0; ipi < Np; ipi++)
                    {
                        for (int istate = 0; istate < num_states; istate++)
                        {
                            state_on_quad(istate) += Phi[ig][ipi] * states_in_element(ipi)(istate);
                        }
                    }
                    ublas::matrix<double> analytical_flux = euler::CalcAnalyticalFlux(state_on_quad, gamma);
                    ublas::vector<double> GPhi_on_quad (2, 0.0);
                    GPhi_on_quad(0) = GPhi[ig][ip][0]; GPhi_on_quad(1) = GPhi[ig][ip][1];
                    ublas::vector<double> GPhi_matmul_invJ (2, 0.0);
                    ublas::vector<double> GPhi_matmul_invJ_matmul_flux (4, 0.0);
                    ublas::axpy_prod(GPhi_on_quad, inv_jacobian, GPhi_matmul_invJ, true);
                    ublas::axpy_prod(GPhi_matmul_invJ, ublas::trans(analytical_flux), GPhi_matmul_invJ_matmul_flux, true);
                    temp_sum += GPhi_matmul_invJ_matmul_flux * det_jacobian * w_quad_2d(ig);
                }
                for (int istate = 0; istate < num_states; istate++)
                {
                    // the contribution from interior, !!! substracted !!!
                    Residual(ielem * Np * num_states + ip * num_states + istate) -= temp_sum(istate);
                }
            }
        }
        else // The curved element
        {
            /*********************************************/
            /* Interior Contribution from Curved element */
            /*********************************************/
            for (int ip = 0; ip < Np; ip++) // loop over all lagrange nodes
            {
                ublas::vector<double> temp_sum (num_states, 0.0); // the temp summation for the quadrature integration
                for (int ig = 0; ig < n_quad_2d; ig++) // loop over all quadrature points
                {
                    // interpolate the state to quadrature points
                    ublas::vector<double> state_on_quad(num_states, 0.0);
                    for (int ipi = 0; ipi < Np; ipi++)
                    {
                        for (int istate = 0; istate < num_states; istate++)
                        {
                            state_on_quad(istate) += Phi[ig][ipi] * states_in_element(ipi)(istate);
                        }
                    }
                    ublas::matrix<double> jacobian = geometry::CalcJacobianCurved(mesh, ielem, GPhi, n_quad_2d, ig);
                    ublas::matrix<double> inv_jacobian (jacobian.size1(), jacobian.size2());
                    double det_jacobian = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
                    InvertMatrix(jacobian, inv_jacobian);

                    ublas::matrix<double> analytical_flux = euler::CalcAnalyticalFlux(state_on_quad, gamma);
                    ublas::vector<double> GPhi_on_quad (2, 0.0);
                    GPhi_on_quad(0) = GPhi[ig][ip][0]; GPhi_on_quad(1) = GPhi[ig][ip][1];
                    ublas::vector<double> GPhi_matmul_invJ (2, 0.0);
                    ublas::vector<double> GPhi_matmul_invJ_matmul_flux (4, 0.0);
                    ublas::axpy_prod(GPhi_on_quad, inv_jacobian, GPhi_matmul_invJ, true);
                    ublas::axpy_prod(GPhi_matmul_invJ, ublas::trans(analytical_flux), GPhi_matmul_invJ_matmul_flux, true);
                    temp_sum += GPhi_matmul_invJ_matmul_flux * det_jacobian * w_quad_2d(ig);
                }
                for (int istate = 0; istate < num_states; istate++)
                {
                    // the contribution from interior, !!! substracted !!!
                    Residual(ielem * Np * num_states + ip * num_states + istate) -= temp_sum(istate);
                }
                if (ielem == 71 && ip == 0)
                    std::cout << "interior:" << ublas::norm_2(temp_sum) << std::endl;
            }
        }
    } // End Loop over elements

    // Loop through interior edges, calculate the edge flux
    for (int iedge = 0; iedge < mesh.I2E.size(); iedge++)
    {
        int ielemL = mesh.I2E[iedge][0] - 1; int ielemR = mesh.I2E[iedge][2] - 1;
        int ilocL = mesh.I2E[iedge][1] - 1; int ilocR = mesh.I2E[iedge][3] - 1;
        ublas::vector<ublas::vector<double> > uL(Np, ublas::vector<double> (num_states, 0.0));
        ublas::vector<ublas::vector<double> > uR(Np, ublas::vector<double> (num_states, 0.0));
        ublas::vector<double> norm_vec(2, 0.0);
        norm_vec(0) = mesh.In[iedge][0];
        norm_vec(1) = mesh.In[iedge][1];
        double jacobian_edge = mesh.In[iedge][2];

        for (int ip = 0; ip < Np; ip++)
        {
            ublas::vector<double> stateL(num_states);
            ublas::vector<double> stateR(num_states);
            for (int istate = 0; istate < num_states; istate++)
            {
                stateL(istate) = States(ielemL * Np * num_states + ip * num_states + istate);
                stateR(istate) = States(ielemR * Np * num_states + ip * num_states + istate);
            }
            uL(ip) = stateL;
            uR(ip) = stateR;
        }
        for (int ip = 0; ip < Np; ip++)
        {
            ublas::vector<double> temp_sum_R (num_states, 0.0);
            ublas::vector<double> temp_sum_L (num_states, 0.0); // the temp summation for the quadrature integration
            // Now do the integration using 1d quad points
            for (int ig = 0; ig < n_quad_1d; ig++)
            {
                // interpolate the LEFT and RIGHT state to a certain quadrature point
                ublas::vector<double> uL_quad (num_states, 0.0), uR_quad (num_states, 0.0);
                for (int ipi = 0; ipi < Np; ipi++)
                {
                    uL_quad += Phi_1D[ilocL][ig][ipi] * uL(ipi);
                    uR_quad += Phi_1D[ilocR][n_quad_1d - 1 - ig][ipi] * uR(ipi); // Note that the
                }
                double mws = 0.0;
                ublas::vector<double> numerical_flux = euler::CalcNumericalFlux(uL_quad, uR_quad, norm_vec, gamma, "roe", mws);
                temp_sum_L += Phi_1D[ilocL][ig][ip] * numerical_flux * jacobian_edge * w_quad_1d(ig);
                temp_sum_R -= Phi_1D[ilocR][n_quad_1d - 1 - ig][ip] * numerical_flux * jacobian_edge * w_quad_1d(ig);
            }
            for (int istate = 0; istate < num_states; istate++)
            {
                // the contribution from edge, !!! ADD !!!
                Residual(ielemL * Np * num_states + ip * num_states + istate) += temp_sum_L(istate);
                Residual(ielemR * Np * num_states + ip * num_states + istate) += temp_sum_R(istate);
            }
            if (ielemL == 71 && ip == 0)
                std::cout << "L " << ublas::norm_2(temp_sum_L) << std::endl;
            if (ielemR == 71 && ip == 0)
                std::cout << "R " << ublas::norm_2(temp_sum_R) << std::endl;
        }
    }

    // Loop through the boundary edges
    for (int iedge = 0; iedge < mesh.B2E.size(); iedge++)
    {
        int ielemL = mesh.B2E[iedge][0] - 1;
        int ilocL = mesh.B2E[iedge][1] - 1;
        ublas::vector<ublas::vector<double> > uL(Np, ublas::vector<double> (num_states, 0.0));
        char* boundary_type;

        // Get the boundary type
        switch (mesh.B2E[iedge][2])
        {
            case 1:
            boundary_type = param.bound0;
                break;
            case 2:
            boundary_type = param.bound1;
                break;
            case 3:
            boundary_type = param.bound2;
                break;
            case 4:
            boundary_type = param.bound3;
                break;
        }
        for (int ip = 0; ip < Np; ip++)
        {
            ublas::vector<double> stateL(num_states);
            ublas::vector<double> stateR(num_states);
            for (int istate = 0; istate < num_states; istate++)
            {
                stateL(istate) = States(ielemL * Np * num_states + ip * num_states + istate);
            }
            uL(ip) = stateL;
        } // Got the inside state vector
        if (mesh.Bn[iedge][3] != 1) // not a curved element
        {
            ublas::vector<double> norm_vec(2, 0.0);
            norm_vec(0) = mesh.Bn[iedge][0];
            norm_vec(1) = mesh.Bn[iedge][1];
            double jacobian_edge = mesh.Bn[iedge][2];
            for (int ip = 0; ip < Np; ip++)
            {
                ublas::vector<double> temp_sum_L (num_states, 0.0); // the temp summation for the quadrature integration
                // Now do the integration using 1d quad points
                for (int ig = 0; ig < n_quad_1d; ig++)
                {
                    // interpolate the LEFT and RIGHT state to a certain quadrature point
                    ublas::vector<double> uL_quad (num_states, 0.0);
                    for (int ipi = 0; ipi < Np; ipi++)
                    {
                        uL_quad += Phi_1D[ilocL][ig][ipi] * uL(ipi);
                    }
                    double mws = 0.0;
                    // Apply the boudary condition
                    ublas::vector<double> numerical_flux = euler::ApplyBoundaryCondition(uL_quad, norm_vec, boundary_type, param, mws);
                    temp_sum_L += Phi_1D[ilocL][ig][ip] * numerical_flux * jacobian_edge * w_quad_1d(ig);
                }
                for (int istate = 0; istate < num_states; istate++)
                {
                    // the contribution from edge, !!! ADD !!!
                    Residual(ielemL * Np * num_states + ip * num_states + istate) += temp_sum_L(istate);
                }
            }
        }
        else // is a curved edge
        {
            // Get the geometry points on the edge
            ublas::vector<ublas::vector<double> > edge_coord = geometry::GetEdgeCoordinates(mesh, iedge);
            ublas::vector<int> edge_coord_ind = geometry::GetEdgeCoordinatesIndex(mesh, iedge);
            ublas::vector<ublas::vector<double> > norm_on_quad_curved(n_quad_1d, ublas::vector<double> (2));
            for (int ig = 0; ig < n_quad_1d; ig++)
            {
                ublas::vector<double> tangent(2, 0.0);
                for (int iq = 0; iq < q + 1; iq++)
                {
                    int local_lagrange_ind = edge_coord_ind(iq);
                    tangent(0) += edge_coord(iq)(0) * GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][0];
                    tangent(1) += edge_coord(iq)(1) * GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][1];
                }
                norm_on_quad_curved(ig)(0) = tangent(1);
                norm_on_quad_curved(ig)(1) = -1.0 * tangent(0);
            }

            for (int ip = 0; ip < Np; ip++)
            {
                ublas::vector<double> temp_sum_L (num_states, 0.0); // the temp summation for the quadrature integration
                // Now do the integration using 1d quad points
                for (int ig = 0; ig < n_quad_1d; ig++)
                {
                    // interpolate the LEFT and RIGHT state to a certain quadrature point
                    ublas::vector<double> uL_quad (num_states, 0.0);
                    for (int ipi = 0; ipi < Np; ipi++)
                    {
                        uL_quad += Phi_1D[ilocL][ig][ipi] * uL(ipi);
                    }
                    double mws = 0.0;
                    // Apply the boudary condition
                    double jacobian_edge = ublas::norm_2(norm_on_quad_curved(ig));
                    ublas::vector<double> norm_vec = norm_on_quad_curved(ig) / jacobian_edge;
                    ublas::vector<double> numerical_flux = euler::ApplyBoundaryCondition(uL_quad, norm_vec, boundary_type, param, mws);
                    temp_sum_L += Phi_1D[ilocL][ig][ip] * numerical_flux * jacobian_edge * w_quad_1d(ig);
                    std::cout << jacobian_edge << std::endl;
                }
                std::cout << "================" << std::endl;
                for (int istate = 0; istate < num_states; istate++)
                {
                    // the contribution from edge, !!! ADD !!!
                    Residual(ielemL * Np * num_states + ip * num_states + istate) += temp_sum_L(istate);
                }
                if (ielemL == 71 && ip == 0)
                    std::cout << "boundary:" << ublas::norm_2(temp_sum_L) << std::endl;
            }

        }

    }
    return Residual;
}


ResData CalcResData(TriMesh mesh, int p)
{
    ResData resdata;

    int Np = int((p + 1) * (p + 2) / 2);
    int Nq = mesh.E[mesh.CurvedIndex[0]].size();
    int q = int((sqrt(1 + 8.0 * Nq) - 3) / 2);
    ublas::matrix<double> TriLagrangeCoeff = lagrange::TriangleLagrange2D(p);
    ublas::matrix<double> TriLagrangeCoeff_Curved = lagrange::TriangleLagrange2D(q);

    int n_quad_1d, n_quad_2d;
    std::vector<double> x_quad_1d_std, w_quad_1d_std, x_quad_2d_std, w_quad_2d_std;
    GetQuadraturePointsWeight1D(9, n_quad_1d,  x_quad_1d_std, w_quad_1d_std); // order 9 is sufficient
    GetQuadraturePointsWeight2D(9, n_quad_2d, x_quad_2d_std, w_quad_2d_std);
    ublas::vector<double> x_quad_1d = utils::StdToBoostVector(x_quad_1d_std);
    ublas::vector<double> w_quad_1d = utils::StdToBoostVector(w_quad_1d_std);
    ublas::vector<double> x_quad_2d = utils::StdToBoostVector(x_quad_2d_std);
    ublas::vector<double> w_quad_2d = utils::StdToBoostVector(w_quad_2d_std);

    resdata.n_quad_1d = n_quad_1d;
    resdata.n_quad_2d = n_quad_2d;
    resdata.x_quad_1d = x_quad_1d;
    resdata.w_quad_1d = w_quad_1d;
    resdata.x_quad_2d = x_quad_2d;
    resdata.w_quad_2d = w_quad_2d;
    resdata.Np = Np;
    resdata.Nq = Nq;

    arr_2d Phi(boost::extents[n_quad_2d][Np]); resdata.Phi.resize(boost::extents[n_quad_2d][Np]);
    arr_3d GPhi(boost::extents[n_quad_2d][Np][2]); resdata.GPhi.resize(boost::extents[n_quad_2d][Np][2]);
    // Pre-calculate the basis functions on edge quad nodes
    arr_3d Phi_1D(boost::extents[3][n_quad_1d][Np]); resdata.Phi_1D.resize(boost::extents[3][n_quad_1d][Np]);
    arr_3d Phi_1D_Curved(boost::extents[3][n_quad_1d][Nq]); resdata.Phi_1D_Curved.resize(boost::extents[3][n_quad_1d][Nq]);
    arr_4d GPhi_1D(boost::extents[3][n_quad_1d][Np][2]); resdata.GPhi_1D.resize(boost::extents[3][n_quad_1d][Np][2]);
    arr_4d GPhi_1D_Curved(boost::extents[3][n_quad_1d][Nq][2]); resdata.GPhi_1D_Curved.resize(boost::extents[3][n_quad_1d][Nq][2]);

    for (int ig = 0; ig < n_quad_2d; ig++) // for each quadrature point
    {
        double xi = x_quad_2d[2 * ig];
        double eta = x_quad_2d[2 * ig + 1];
        ublas::matrix<double> gphi = lagrange::CalcBaseFunctionGradient(TriLagrangeCoeff, xi, eta);
        ublas::vector<double> phi = lagrange::CalcBaseFunction(TriLagrangeCoeff, xi, eta);
        for (int ipi = 0; ipi < Np; ipi++)
        {
            Phi[ig][ipi] = phi(ipi);
            GPhi[ig][ipi][0] = gphi(ipi, 0);
            GPhi[ig][ipi][1] = gphi(ipi, 1);
        }
    }

    for (int num_edge = 0; num_edge < 3; num_edge++)
    {
        for (int ig = 0; ig < n_quad_1d; ig++) // for each quadrature point
        {
            double xi, eta;
            switch (num_edge)
            {
                case 0:
                    xi = 1 - x_quad_1d[ig]; eta = 1 - xi;
                break;
                case 1:
                    xi = 0.0;           eta = 1- x_quad_1d[ig];
                break;
                case 2:
                    xi = x_quad_1d[ig]; eta = 0.0;
                break;
            }
            ublas::vector<double> phi = lagrange::CalcBaseFunction(TriLagrangeCoeff, xi, eta);
            ublas::vector<double> phi_curved = lagrange::CalcBaseFunction(TriLagrangeCoeff_Curved, xi, eta);
            ublas::matrix<double> gphi = lagrange::CalcBaseFunctionGradient(TriLagrangeCoeff, xi, eta);
            ublas::matrix<double> gphi_curved = lagrange::CalcBaseFunctionGradient(TriLagrangeCoeff_Curved, xi, eta);
            for (int ip = 0; ip < Np; ip++)
            {
                Phi_1D[num_edge][ig][ip] = phi(ip);
                GPhi_1D[num_edge][ig][ip][0] = gphi(ip, 0);
                GPhi_1D[num_edge][ig][ip][1] = gphi(ip, 1);
            }
            for (int iq = 0; iq < Nq; iq++)
            {
                Phi_1D_Curved[num_edge][ig][iq] = phi_curved(iq);
                GPhi_1D_Curved[num_edge][ig][iq][0] = gphi_curved(iq, 0);
                GPhi_1D_Curved[num_edge][ig][iq][1] = gphi_curved(iq, 1);
            }
        }
    }
    resdata.Phi = Phi;
    resdata.GPhi = GPhi;
    resdata.Phi_1D = Phi_1D;
    resdata.Phi_1D_Curved = Phi_1D_Curved;
    resdata.GPhi_1D = GPhi_1D;
    resdata.GPhi_1D_Curved = GPhi_1D_Curved;
    return resdata;
}