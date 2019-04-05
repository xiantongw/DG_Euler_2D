#include "../include/solver.h"

namespace ublas = boost::numeric::ublas;

namespace solver{

    ublas::vector<double> CalcResidual(TriMesh mesh, Param& param, ResData& resdata, ublas::vector<double> States, ublas::vector<double>& dt, int p)
    {
        // Unroll the mesh information
        int num_element = mesh.num_element;
        int num_states = 4; // Four states for 2D Euler equation
        // The number of lagrange nodes in each element
        int Np = int((p + 1) * (p + 2) / 2);
        int Nq = mesh.E[mesh.CurvedElementIndex[0]].size();
        int q = int((sqrt(1 + 8.0 * Nq) - 3) / 2);
        double gamma = param.gamma;
        // The residual vector
        ublas::vector<double> Residual (num_element * Np * num_states, 0.0);
        ublas::vector<double> mws_tally(num_element, 0.0);
        // Pull out the ResData
        arr_2d Phi = resdata.Phi;
        arr_3d GPhi = resdata.GPhi;
        arr_3d GPhi_Curved = resdata.GPhi_Curved;
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

        Residual.clear();
        mws_tally.clear();

        // Loop over elements
        // 1. Construct Jacobian 2.evealute interior contribution to the residual
        // Notice that there are two sets of elements in the mesh, curved and not curved
        for (int ielem_linear = 0; ielem_linear < mesh.LinearElementIndex.size(); ielem_linear++)
        {
            int ielem = mesh.LinearElementIndex[ielem_linear];
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
            /*********************************************/
            /* Interior Contribution from linear element */
            /*********************************************/
            // Construct the Jacobian matrix
            ublas::matrix<double> jacobian = resdata.jacobian_in_linear_elements(ielem);
            double det_jacobian = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
            ublas::matrix<double> inv_jacobian = resdata.invjacobian_in_linear_elements(ielem);
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
                    GPhi_on_quad(0) = GPhi[ig][ip][0];
                    GPhi_on_quad(1) = GPhi[ig][ip][1];
                    ublas::vector<double> GPhi_matmul_invJ (2, 0.0);
                    ublas::vector<double> GPhi_matmul_invJ_matmul_flux (4, 0.0);
                    ublas::axpy_prod(GPhi_on_quad, inv_jacobian, GPhi_matmul_invJ, false);
                    ublas::axpy_prod(GPhi_matmul_invJ, ublas::trans(analytical_flux), GPhi_matmul_invJ_matmul_flux, false);
                    temp_sum += GPhi_matmul_invJ_matmul_flux * det_jacobian * w_quad_2d(ig);
                }
                for (int istate = 0; istate < num_states; istate++)
                {
                    // the contribution from interior, !!! substracted !!!
                    Residual(ielem * Np * num_states + ip * num_states + istate) -= temp_sum(istate);
                }
            }
        }


        for (int ielem_curved = 0; ielem_curved < mesh.CurvedElementIndex.size(); ielem_curved++)
        {
            int ielem = mesh.CurvedElementIndex[ielem_curved];
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
                    ublas::matrix<double> jacobian = resdata.jacobian_in_curved_elements(ielem, ig);
                    double det_jacobian = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
                    ublas::matrix<double> inv_jacobian = resdata.invjacobian_in_curved_elements(ielem, ig);

                    ublas::matrix<double> analytical_flux = euler::CalcAnalyticalFlux(state_on_quad, gamma);
                    ublas::vector<double> GPhi_on_quad (2, 0.0);
                    GPhi_on_quad(0) = GPhi[ig][ip][0]; GPhi_on_quad(1) = GPhi[ig][ip][1];
                    ublas::vector<double> GPhi_matmul_invJ (2, 0.0);
                    ublas::vector<double> GPhi_matmul_invJ_matmul_flux (4, 0.0);
                    ublas::axpy_prod(GPhi_on_quad, inv_jacobian, GPhi_matmul_invJ, false);
                    ublas::axpy_prod(GPhi_matmul_invJ, ublas::trans(analytical_flux), GPhi_matmul_invJ_matmul_flux, false);
                    temp_sum += GPhi_matmul_invJ_matmul_flux * det_jacobian * w_quad_2d(ig);
                }
                for (int istate = 0; istate < num_states; istate++)
                {
                    // the contribution from interior, !!! substracted !!!
                    Residual(ielem * Np * num_states + ip * num_states + istate) -= temp_sum(istate);
                }
            }
        } // End Loop over elements

        // Loop through interior edges, calculate the edge flux
        for (int iedge = 0; iedge < mesh.I2E.size(); iedge++)
        {
            double mws_recorded = 0.0;
            int ielemL = mesh.I2E[iedge][0] - 1; int ielemR = mesh.I2E[iedge][2] - 1;
            int ilocL = mesh.I2E[iedge][1] - 1; int ilocR = mesh.I2E[iedge][3] - 1;
            ublas::vector<ublas::vector<double> > uL(Np, ublas::vector<double> (num_states, 0.0));
            ublas::vector<ublas::vector<double> > uR(Np, ublas::vector<double> (num_states, 0.0));
            ublas::vector<double> norm_vec(2, 0.0);
            norm_vec(0) = mesh.In[iedge][0];
            norm_vec(1) = mesh.In[iedge][1];
            double jacobian_edge = mesh.In[iedge][2];
            double mws = 0.0;
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
                        uR_quad += Phi_1D[ilocR][n_quad_1d - 1 - ig][ipi] * uR(ipi);
                    }
                    ublas::vector<double> numerical_flux = euler::CalcNumericalFlux(uL_quad, uR_quad, norm_vec, gamma, "roe", mws);
                    if (mws_recorded < mws)
                        mws_recorded = mws;
                    temp_sum_L += Phi_1D[ilocL][ig][ip] * numerical_flux * jacobian_edge * w_quad_1d(ig);
                    temp_sum_R -= Phi_1D[ilocR][n_quad_1d - 1 - ig][ip] * numerical_flux * jacobian_edge * w_quad_1d(ig);
                }
                for (int istate = 0; istate < num_states; istate++)
                {
                    // the contribution from edge, !!! ADD !!!
                    Residual(ielemL * Np * num_states + ip * num_states + istate) += temp_sum_L(istate);
                    Residual(ielemR * Np * num_states + ip * num_states + istate) += temp_sum_R(istate);
                }
            }
            mws_tally(ielemL) += mws_recorded * jacobian_edge;
            mws_tally(ielemR) += mws_recorded * jacobian_edge;
        }

        // Loop through the boundary curved edges
        for (int iedge_curved = 0; iedge_curved < mesh.CurvedEdgeIndex.size(); iedge_curved++)
        {
            int iedge = mesh.CurvedEdgeIndex[iedge_curved];
            int ielemL = mesh.B2E[iedge][0] - 1;
            int ilocL = mesh.B2E[iedge][1] - 1;
            ublas::vector<ublas::vector<double> > uL(Np, ublas::vector<double> (num_states, 0.0));
            string boundary_type;
            boundary_type = param.bound0;
            double mws = 0.0, mws_recorded = 0.0;
            double jacobian_edge, jacobian_edge_recorded = 0.0;
            // Get the boundary type
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
                    double deriv_along_edge = 0.0;
                    switch (ilocL)
                    {
                        case 0:
                            tangent(0) += - edge_coord(iq)(0) * GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][0]
                                        + edge_coord(iq)(0) * GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][1];
                            tangent(1) += - edge_coord(iq)(1) * GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][0]
                                        + edge_coord(iq)(1) * GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][1];
                            break;
                        case 1:
                            deriv_along_edge = -GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][1];
                            tangent(0) += edge_coord(iq)(0) * deriv_along_edge;
                            tangent(1) += edge_coord(iq)(1) * deriv_along_edge;
                            break;
                        case 2:
                            deriv_along_edge = GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][0];
                            tangent(0) += edge_coord(iq)(0) * deriv_along_edge;
                            tangent(1) += edge_coord(iq)(1) * deriv_along_edge;
                            break;
                        default:
                            break;
                    }
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
                    // Apply the boudary condition
                    jacobian_edge = ublas::norm_2(norm_on_quad_curved(ig));
                    if (jacobian_edge > jacobian_edge_recorded)
                        jacobian_edge_recorded = jacobian_edge;
                    ublas::vector<double> norm_vec = norm_on_quad_curved(ig) / jacobian_edge;
                    ublas::vector<double> numerical_flux = euler::ApplyBoundaryCondition(uL_quad, norm_vec, boundary_type, param, mws);
                    if (mws_recorded < mws)
                        mws_recorded = mws;
                    temp_sum_L += Phi_1D[ilocL][ig][ip] * numerical_flux * jacobian_edge * w_quad_1d(ig);
                }
                for (int istate = 0; istate < num_states; istate++)
                {
                    // the contribution from edge, !!! ADD !!!
                    Residual(ielemL * Np * num_states + ip * num_states + istate) += temp_sum_L(istate);
                }
            }
            mws_tally(ielemL) += mws_recorded * jacobian_edge_recorded;
        }

        for (int iedge_linear = 0; iedge_linear < mesh.LinearEdgeIndex.size(); iedge_linear++)
        {
            int iedge = mesh.LinearEdgeIndex[iedge_linear];
            int ielemL = mesh.B2E[iedge][0] - 1;
            int ilocL = mesh.B2E[iedge][1] - 1;
            ublas::vector<ublas::vector<double> > uL(Np, ublas::vector<double> (num_states, 0.0));
            string boundary_type;
            double mws = 0.0, mws_recorded = 0.0;
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
                    // Apply the boudary condition
                    ublas::vector<double> numerical_flux = euler::ApplyBoundaryCondition(uL_quad, norm_vec, boundary_type, param, mws);
                    if (mws_recorded < mws)
                        mws_recorded = mws;
                    temp_sum_L += Phi_1D[ilocL][ig][ip] * numerical_flux * jacobian_edge * w_quad_1d(ig);
                }
                for (int istate = 0; istate < num_states; istate++)
                {
                    // the contribution from edge, !!! ADD !!!
                    Residual(ielemL * Np * num_states + ip * num_states + istate) += temp_sum_L(istate);
                }
            }
            mws_tally(ielemL) += mws_recorded * jacobian_edge;
        }
        // Calculate dtA
        for (int i = 0; i < num_element; i++)
        {
            dt(i) = 2.0 * mesh.Area[i] * param.cfl / mws_tally(i);
        }
        return Residual;
    }


    void CalcResData(TriMesh mesh, int p, ResData& resdata)
    {
        int Np = int((p + 1) * (p + 2) / 2);
        int Nq = mesh.E[mesh.CurvedElementIndex[0]].size();
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
        arr_2d Phi_Curved(boost::extents[n_quad_2d][Nq]); resdata.Phi_Curved.resize(boost::extents[n_quad_2d][Nq]);
        arr_3d GPhi(boost::extents[n_quad_2d][Np][2]); resdata.GPhi.resize(boost::extents[n_quad_2d][Np][2]);
        arr_3d GPhi_Curved(boost::extents[n_quad_2d][Nq][2]); resdata.GPhi_Curved.resize(boost::extents[n_quad_2d][Nq][2]);
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

        for (int ig = 0; ig < n_quad_2d; ig++) // for each quadrature point
        {
            double xi = x_quad_2d[2 * ig];
            double eta = x_quad_2d[2 * ig + 1];
            ublas::vector<double> phi_curved = lagrange::CalcBaseFunction(TriLagrangeCoeff_Curved, xi, eta);
            ublas::matrix<double> gphi_curved = lagrange::CalcBaseFunctionGradient(TriLagrangeCoeff_Curved, xi, eta);
            for (int ipi = 0; ipi < Nq; ipi++)
            {
                Phi_Curved[ig][ipi] = phi_curved(ipi);
                GPhi_Curved[ig][ipi][0] = gphi_curved(ipi, 0);
                GPhi_Curved[ig][ipi][1] = gphi_curved(ipi, 1);
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
                        xi = 0.0;           eta = 1 - x_quad_1d[ig];
                    break;
                    case 2:
                        xi = x_quad_1d[ig]; eta = 0.0;
                    break;
                }
                ublas::vector<double> phi = lagrange::CalcBaseFunction(TriLagrangeCoeff, xi, eta);
                ublas::matrix<double> gphi = lagrange::CalcBaseFunctionGradient(TriLagrangeCoeff, xi, eta);
                ublas::vector<double> phi_curved = lagrange::CalcBaseFunction(TriLagrangeCoeff_Curved, xi, eta);
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
        // Calculate the jacobian in curved elements on quadrature nodes
        ublas::matrix<ublas::matrix<double> > jacobian_curved(mesh.E.size(), n_quad_2d, ublas::matrix<double> (2, 2, 0.0));
        ublas::matrix<ublas::matrix<double> > inv_jacobian_curved(mesh.E.size(), n_quad_2d, ublas::matrix<double> (2, 2, 0.0));
        ublas::vector<ublas::matrix<double> > jacobian_linear(mesh.E.size(), ublas::matrix<double> (2, 2, 0.0));
        ublas::vector<ublas::matrix<double> > inv_jacobian_linear(mesh.E.size(), ublas::matrix<double> (2, 2, 0.0));
        for (int ielem = 0; ielem < mesh.E.size(); ielem++)
        {
            if (mesh.isCurved[ielem])
            {
                for (int ig = 0; ig < n_quad_2d; ig++)
                {
                    ublas::matrix<double> jacobian(2, 2, 0.0), inv_jacobian(2, 2, 0.0);
                    jacobian = geometry::CalcJacobianCurved(mesh, ielem, GPhi_Curved, n_quad_2d, ig);
                    inv_jacobian = utils::Invert22Matrix(jacobian);
                    jacobian_curved(ielem, ig) = jacobian;
                    inv_jacobian_curved(ielem, ig) = inv_jacobian;

                }
            }
            else
            {
                ublas::matrix<double> jacobian(2, 2, 0.0), inv_jacobian(2, 2, 0.0);
                jacobian = geometry::CalcJacobianLinear(mesh, ielem);
                inv_jacobian = utils::Invert22Matrix(jacobian);
                jacobian_linear(ielem) = jacobian;
                inv_jacobian_linear(ielem) = inv_jacobian;
            }
        }

        resdata.jacobian_in_curved_elements = jacobian_curved;
        resdata.jacobian_in_linear_elements = jacobian_linear;
        resdata.invjacobian_in_curved_elements = inv_jacobian_curved;
        resdata.invjacobian_in_linear_elements = inv_jacobian_linear;
        resdata.Phi = Phi;
        resdata.Phi_Curved = Phi_Curved;
        resdata.GPhi = GPhi;
        resdata.GPhi_Curved = GPhi_Curved;
        resdata.Phi_1D = Phi_1D;
        resdata.Phi_1D_Curved = Phi_1D_Curved;
        resdata.GPhi_1D = GPhi_1D;
        resdata.GPhi_1D_Curved = GPhi_1D_Curved;

    }

        ublas::vector<double> TimeMarching_TVDRK3(TriMesh mesh, Param& param, ResData& resdata, ublas::vector<double> States_old, ublas::vector<ublas::matrix<double> > invM, int p, int& converged, double& norm_residual)
    {
        ublas::vector<double> States_new = States_old;
        ublas::vector<double> States_1 = States_old * 0.0;
        ublas::vector<double> States_2 = States_old * 0.0;
        ublas::vector<double> Residual_1(States_old.size(), 0.0);
        ublas::vector<double> Residual_2(States_old.size(), 0.0);
        int num_elements = invM.size(); int num_states = 4;
        int Np = int((p + 1) * (p + 2) / 2);
        double eps = param.eps;
        converged = 0;
        ublas::vector<double> dt (num_elements, 0.0), dt_temp (num_elements, 0.0);
        ublas::vector<double> Residual = CalcResidual(mesh, param, resdata, States_old, dt, p); // Caculate the residual, and the time step
        // Caculate the 1st state in TVDRK3, the first step
        for (int ielem = 0; ielem < num_elements; ielem++)
        {
            // Get the states in this element
            ublas::matrix<double> u(Np, num_states, 0.0);
            ublas::matrix<double> R(Np, num_states, 0.0);
            ublas::matrix<double> u_1(Np, num_states, 0.0);
            for (int ip = 0; ip < Np; ip++)
            {
                for (int istate = 0; istate < num_states; istate++)
                {
                    u(ip, istate) = States_old(ielem * Np * num_states + ip * num_states + istate);
                    R(ip, istate) = Residual(ielem * Np * num_states + ip * num_states + istate);
                }
            }
            ublas::matrix<double> invM_mul_R(Np, num_states, 0.0);
            ublas::axpy_prod(invM(ielem), R, invM_mul_R, false);
            u_1 = u - dt(ielem) * invM_mul_R;
            // Apply the values of u_FE to the long vector State_FE
            for (int ip = 0; ip < Np; ip++)
            {
                for (int istate = 0; istate < num_states; istate++)
                {
                    States_1(ielem * Np * num_states + ip * num_states + istate) = u_1(ip, istate);
                }
            }
        }
        Residual_1  = CalcResidual(mesh, param, resdata, States_1, dt_temp, p);
        // The second step of RK3
        for (int ielem = 0; ielem < num_elements; ielem++)
        {
            // check if an element is converged
            // Get the states in this element
            ublas::matrix<double> u(Np, num_states, 0.0);
            ublas::matrix<double> R(Np, num_states, 0.0);
            ublas::matrix<double> u_1(Np, num_states, 0.0);
            ublas::matrix<double> u_2(Np, num_states, 0.0);
            for (int ip = 0; ip < Np; ip++)
            {
                for (int istate = 0; istate < num_states; istate++)
                {
                    u(ip, istate) = States_old(ielem * Np * num_states + ip * num_states + istate);
                    u_1(ip, istate) = States_1(ielem * Np * num_states + ip * num_states + istate);
                    R(ip, istate) = Residual_1(ielem * Np * num_states + ip * num_states + istate);
                }
            }
            ublas::matrix<double> invM_mul_R(Np, num_states, 0.0);
            ublas::axpy_prod(invM(ielem), R, invM_mul_R, false);
            u_2 = 0.75 * u + 0.25 * u_1  - 0.25 * dt(ielem) * invM_mul_R;
            // Apply the values of u_FE to the long vector State_FE
            for (int ip = 0; ip < Np; ip++)
            {
                for (int istate = 0; istate < num_states; istate++)
                {
                    States_2(ielem * Np * num_states + ip * num_states + istate) = u_2(ip, istate);
                }
            }
        }
        Residual_2  = CalcResidual(mesh, param, resdata, States_2, dt_temp, p);
        for (int ielem = 0; ielem < num_elements; ielem++)
        {
            // check if an element is converged
            // Get the states in this element
            ublas::matrix<double> u(Np, num_states, 0.0);
            ublas::matrix<double> R(Np, num_states, 0.0);
            ublas::matrix<double> u_1(Np, num_states, 0.0);
            ublas::matrix<double> u_2(Np, num_states, 0.0);
            ublas::matrix<double> u_new(Np, num_states, 0.0);
            for (int ip = 0; ip < Np; ip++)
            {
                for (int istate = 0; istate < num_states; istate++)
                {
                    u(ip, istate) = States_old(ielem * Np * num_states + ip * num_states + istate);
                    u_1(ip, istate) = States_1(ielem * Np * num_states + ip * num_states + istate);
                    u_2(ip, istate) = States_2(ielem * Np * num_states + ip * num_states + istate);
                    u_new(ip, istate) = States_new(ielem * Np * num_states + ip * num_states + istate);
                    R(ip, istate) = Residual_2(ielem * Np * num_states + ip * num_states + istate);
                }
            }
            ublas::matrix<double> invM_mul_R(Np, num_states, 0.0);
            ublas::axpy_prod(invM(ielem), R, invM_mul_R, false);
            u_new = (1.0 / 3) * u + (2.0 / 3) * u_2 - (2.0 / 3) * dt(ielem) * invM_mul_R;
            // Apply the values of u_new to the long vector State_new
            for (int ip = 0; ip < Np; ip++)
            {
                for (int istate = 0; istate < num_states; istate++)
                {
                    States_new(ielem * Np * num_states + ip * num_states + istate) = u_new(ip, istate);
                }
            }
        }
        ublas::vector<double> Residual_new  = CalcResidual(mesh, param, resdata, States_new, dt_temp, p);
        norm_residual = ublas::norm_inf(Residual_new);
        if (norm_residual < eps)
            converged = 1;
        return  States_new;
    }

    void PostProc(TriMesh mesh, ublas::vector<double> States, int p, ublas::vector<ublas::matrix<double> >& Nodes, ublas::vector<ublas::matrix<double> >& States_on_Nodes)
    {
        int num_elements = mesh.E.size();
        int Np = int((p + 1) * (p + 2) / 2);
        int num_states = 4;
        ublas::vector<ublas::vector<double> > node_physical;
        if (p == 0)
        {
            int Np_solution = 3;
            for (int ielem = 0; ielem < num_elements; ielem++)
            {
                node_physical = lagrange::MapReferenceToPhysical(mesh, ielem, 1, geometry::BumpFunction);
                for (int ip = 0; ip < Np_solution; ip++)
                {
                    for (int i = 0; i < 2; i++)
                    {
                        Nodes(ielem)(ip, i) = node_physical(ip)(i);
                    }
                    for (int istate = 0; istate < num_states; istate++)
                    {
                        int ipi = int (ip / 3);
                        States_on_Nodes(ielem)(ip, istate) = States(ielem * Np * num_states + ipi * num_states + istate);
                    }
                }
            }
        }
        else
        {
            for (int ielem = 0; ielem < num_elements; ielem++)
            {
                node_physical = lagrange::MapReferenceToPhysical(mesh, ielem, p, geometry::BumpFunction);
                for (int ip = 0; ip < Np; ip++)
                {
                    for (int i = 0; i < 2; i++)
                    {
                        Nodes(ielem)(ip, i) = node_physical(ip)(i);
                    }
                    for (int istate = 0; istate < num_states; istate++)
                    {
                        States_on_Nodes(ielem)(ip, istate) = States(ielem * Np * num_states + ip * num_states + istate);
                    }
                }
            }
        }
    }

    void CalcScalarOutputs(TriMesh mesh, ResData resdata, ublas::vector<ublas::matrix<double> >& States_on_Nodes, ublas::vector<ublas::matrix<double> >& Nodes,
                            Param param, double& err_entropy, double& coeff_lift, double& coeff_drag, std::vector<std::vector<double> >& p_coeff_dist)
    {
        // The constants
        double p_inf = param.p_inf;
        double m_inf = param.mach_inf;
        double R = param.R;
        double gamma = param.gamma;
        double h = param.h;
        // Calculate the derived quantities from the free-stream state
        double T_t = 1.0 + 0.5 * (gamma - 1.0) * m_inf * m_inf;
        double p_t = pow(T_t, gamma / (gamma - 1.0));
        double rho_t = p_t / (R * T_t);
        double s_t = p_t / pow(rho_t, gamma);
        // Caculate pressure and stagnation entropy on nodes
        int num_element = mesh.E.size();
        int Np = States_on_Nodes(0).size1();
        int p = int((sqrt(1 + 8.0 * Np) - 3) / 2);
        ublas::vector<ublas::matrix<double> > p_on_nodes(num_element, ublas::matrix<double>(Np, 1, 0.0));
        ublas::vector<ublas::matrix<double> > rel_s2_on_nodes(num_element, ublas::matrix<double>(Np, 1, 0.0));
        for (int ielem = 0; ielem < num_element; ielem++)
        {
            for (int ip = 0; ip < Np; ip++)
            {
                double rho = States_on_Nodes(ielem)(ip, 0);
                double rhou = States_on_Nodes(ielem)(ip, 1);
                double rhov = States_on_Nodes(ielem)(ip, 2);
                double rhoE = States_on_Nodes(ielem)(ip, 3);
                p_on_nodes(ielem)(ip, 0) = (gamma - 1.0) * (rhoE - 0.5 * (rhou * rhou + rhov * rhov) / rho);
                double s = p_on_nodes(ielem)(ip, 0) / pow(rho, gamma);
                rel_s2_on_nodes(ielem)(ip, 0) = (s / s_t - 1) * (s / s_t - 1);
            }
        }
        double total_area = 0.0;
        err_entropy = 0.0;
        // Calculate the entropy error
        // linear elements
        for (int i_linear_elem = 0; i_linear_elem < mesh.LinearElementIndex.size(); i_linear_elem++)
        {
            int ielem = mesh.LinearElementIndex[i_linear_elem];
            ublas::matrix<double> jacobian = resdata.jacobian_in_linear_elements(ielem);
            double det_jacobian = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
            total_area += det_jacobian;
            for (int ig = 0; ig < resdata.n_quad_2d; ig++) // integrate on quadrature points
            {
                for (int ip = 0; ip < Np; ip++)
                {
                    err_entropy += resdata.Phi[ig][ip] * rel_s2_on_nodes(ielem)(ip, 0) * resdata.w_quad_2d(ig) * det_jacobian;
                }
            }
        }
        // curved elements
        for (int i_curved_elem = 0; i_curved_elem < mesh.CurvedElementIndex.size(); i_curved_elem++)
        {
            int ielem = mesh.CurvedElementIndex[i_curved_elem];
            for (int ig = 0; ig < resdata.n_quad_2d; ig++)
            {
                ublas::matrix<double> jacobian = resdata.jacobian_in_curved_elements(ielem, ig);
                double det_jacobian = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
                total_area += resdata.w_quad_2d(ig) * det_jacobian;
                for (int ip = 0; ip < Np; ip++)
                {
                    err_entropy += resdata.Phi[ig][ip] * rel_s2_on_nodes(ielem)(ip, 0) * resdata.w_quad_2d(ig) * det_jacobian;
                }
            }
        }
        err_entropy = sqrt(err_entropy / total_area); // Entropy error calculated
        // Caculate drag and lift coefficients
        // Loop through the boundary curved edges

        coeff_lift = 0.0, coeff_drag = 0.0;
        for (int iedge_curved = 0; iedge_curved < mesh.CurvedEdgeIndex.size(); iedge_curved++)
        {
            int iedge = mesh.CurvedEdgeIndex[iedge_curved];
            int ielemL = mesh.B2E[iedge][0] - 1;
            int ilocL = mesh.B2E[iedge][1] - 1;
            // Get the geometry points on the edge
            ublas::vector<ublas::vector<double> > edge_coord = geometry::GetEdgeCoordinates(mesh, iedge);
            ublas::vector<int> edge_coord_ind = geometry::GetEdgeCoordinatesIndex(mesh, iedge);
            ublas::vector<ublas::vector<double> > norm_on_quad_curved(resdata.n_quad_1d, ublas::vector<double> (2));
            double Nq = resdata.Nq;
            int q = int((sqrt(1 + 8.0 * Nq) - 3) / 2);
            ublas::vector<double> quad_x_total(resdata.n_quad_1d, 0.0);
            for (int ig = 0; ig < resdata.n_quad_1d; ig++)
            {
                ublas::vector<double> tangent(2, 0.0);
                double quad_x = 0.0;
                for (int iq = 0; iq < q + 1; iq++)
                {
                    int local_lagrange_ind = edge_coord_ind(iq);
                    double deriv_along_edge = 0.0;
                    switch (ilocL)
                    {
                        case 0:
                            tangent(0) += - edge_coord(iq)(0) * resdata.GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][0]
                                        + edge_coord(iq)(0) * resdata.GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][1];
                            tangent(1) += - edge_coord(iq)(1) * resdata.GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][0]
                                        + edge_coord(iq)(1) * resdata.GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][1];
                            quad_x += - edge_coord(iq)(0) * resdata.Phi_1D_Curved[ilocL][ig][local_lagrange_ind];
                            break;
                        case 1:
                            deriv_along_edge = - resdata.GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][1];
                            tangent(0) += edge_coord(iq)(0) * deriv_along_edge;
                            tangent(1) += edge_coord(iq)(1) * deriv_along_edge;
                            quad_x += - edge_coord(iq)(0) * resdata.Phi_1D_Curved[ilocL][ig][local_lagrange_ind];
                            break;
                        case 2:
                            deriv_along_edge = resdata.GPhi_1D_Curved[ilocL][ig][local_lagrange_ind][0];
                            tangent(0) += edge_coord(iq)(0) * deriv_along_edge;
                            tangent(1) += edge_coord(iq)(1) * deriv_along_edge;
                            quad_x += - edge_coord(iq)(0) * resdata.Phi_1D_Curved[ilocL][ig][local_lagrange_ind];
                            break;
                        default:
                            break;
                    }
                }
                norm_on_quad_curved(ig)(0) = tangent(1);
                norm_on_quad_curved(ig)(1) = -1.0 * tangent(0);
                quad_x_total(ig) = quad_x;
            }

            // Now do the integration using 1d quad points
            for (int ig = 0; ig < resdata.n_quad_1d; ig++)
            {
                // interpolate the pressure to quadrature nodes
                double p_quad = 0.0;
                for (int ipi = 0; ipi < Np; ipi++)
                {
                    p_quad += resdata.Phi_1D[ilocL][ig][ipi] * p_on_nodes(ielemL)(ipi, 0);
                }
                double jacobian_edge = ublas::norm_2(norm_on_quad_curved(ig));
                ublas::vector<double> norm_vec = norm_on_quad_curved(ig) / jacobian_edge;
                coeff_lift += (p_quad - p_inf) * norm_vec(1) * jacobian_edge * resdata.w_quad_1d(ig);
                coeff_drag += (p_quad - p_inf) * norm_vec(0) * jacobian_edge * resdata.w_quad_1d(ig);
                double p_coeff = (p_quad - p_inf) / (0.5 * gamma * p_inf * m_inf * m_inf);
                std::vector<double> p_coeff_on_node = {quad_x_total(ig), p_coeff};
                p_coeff_dist.push_back(p_coeff_on_node);
            }
        }
        coeff_lift = coeff_lift / (0.5 * gamma * p_inf * m_inf * m_inf * h);
        coeff_drag = coeff_drag / (0.5 * gamma * p_inf * m_inf * m_inf * h);
    }

} // end namespace solver
