#include "../include/euler.h"

namespace ublas = boost::numeric::ublas;

namespace euler{

    ublas::matrix<double> CalcAnalyticalFlux(ublas::vector<double> state, double gamma)
    {
        /*
            The function caculating the analytical flux from the state vector
            for 2D euler equations
        */
        int nstate = state.size();
        ublas::matrix<double> F(nstate, 2);
        double rho = state(0);
        double u = state(1) / rho;
        double v = state(2) / rho;
        double E = state(3) / rho;
        double q = sqrt(u * u + v * v);
        double p = (gamma - 1) * (rho * E - 0.5 * rho * q * q);
        double H = E + p / rho;
        // assign values to F
        F(0, 0) = rho * u;         F(0, 1) = rho * v;
        F(1, 0) = rho * u * u + p; F(1, 1) = rho * v * u;
        F(2, 0) = rho * u * v;     F(2, 1) = rho * v * v + p;
        F(3, 0) = rho * u * H;     F(3, 1) = rho * v * H;
        return F;
    }

    ublas::vector<double> CalcNumericalFlux(ublas::vector<double> uL, ublas::vector<double> uR, ublas::vector<double> n,
                                            double gamma, char* type_flux, double& mws)
    {
        if (strcasecmp(type_flux, "roe") == 0)
        {
            ublas::vector<double> F_hat(4, 0.0);

            double v_roe_vec[2]; double u_roe;
            double H_roe, c_roe, q_roe;
            double pL, pR, HL, HR, cL, cR, rhoL, rhoR, qL, qR;
            double lambda_1, lambda_2, lambda_3;
            double s1, s2, G1, G2, C1, C2;

            double drho, drhoE, drhov_vec[2];

            //Primitive variables
            rhoL = uL[0];
            rhoR = uR[0];
            pL = (gamma - 1) * (uL[3] - 0.5 * (uL[1]*uL[1] + uL[2]*uL[2]) / uL[0]);
            HL = uL[3] / uL[0] + pL / uL[0];
            pR = (gamma - 1) * (uR[3] - 0.5 * (uR[1]*uR[1] + uR[2]*uR[2]) / uR[0]);

            HR = uR[3] / uR[0] + pR / uR[0];
            cL = sqrt(gamma * pL / rhoL);
            cR = sqrt(gamma * pR / rhoR);
            qL = sqrt(pow(uL[1] / uL[0], 2) + pow(uL[2] / uL[0], 2));
            qR = sqrt(pow(uR[1] / uR[0], 2) + pow(uR[2] / uR[0], 2));
            // calculate v_roe
            v_roe_vec[0] = (sqrt(rhoL) * (uL[1] / uL[0]) + sqrt(rhoR) * (uR[1] / uR[0])) / (sqrt(rhoL) + sqrt(rhoR));
            v_roe_vec[1] = (sqrt(rhoL) * (uL[2] / uL[0]) + sqrt(rhoR) * (uR[2] / uR[0])) / (sqrt(rhoL) + sqrt(rhoR));
            u_roe = v_roe_vec[0] * n[0] + v_roe_vec[1] * n[1];
            // calculate H_roe
            H_roe = (sqrt(rhoL) * HL + sqrt(rhoR) * HR) / (sqrt(rhoL) + sqrt(rhoR));

            // velocity magnitude
            q_roe = sqrt(v_roe_vec[0] * v_roe_vec[0] + v_roe_vec[1] * v_roe_vec[1]);
            // sound speed
            c_roe = sqrt((gamma - 1.0) * (H_roe - 0.5 * (v_roe_vec[0] * v_roe_vec[0] + v_roe_vec[1] * v_roe_vec[1])));
            // eigenvalues
            double lambda[3];
            lambda[0] = u_roe + c_roe;
            lambda[1] = u_roe - c_roe;
            lambda[2] = u_roe;

            // entropy fix for lambda_1, lambda_2 and lambda_3
            double epsilon = 0.1 * c_roe;
            int i;
            for (i = 0; i < 3; i++){
                if (lambda[i] < epsilon && lambda[i] > -epsilon){
                lambda[i] = 0.5 * (epsilon + lambda[i] * lambda[i] / epsilon);
                }
            }
            lambda_1 = lambda[0]; lambda_2 = lambda[1]; lambda_3 = lambda[2];

            // calculate the delta quantities for roe flux
            drho = uR[0] - uL[0];
            drhoE = uR[3] - uL[3];
            drhov_vec[0] = uR[1] - uL[1];
            drhov_vec[1] = uR[2] - uL[2];
            // temporary vars for calculating the Roe Flux
            s1 = 0.5*(fabs(lambda_1) + fabs(lambda_2));
            s2 = 0.5*(fabs(lambda_1) - fabs(lambda_2));
            G1 = (gamma - 1.0) * (q_roe * q_roe * drho / 2 + drhoE - (v_roe_vec[0] * drhov_vec[0] + v_roe_vec[1] * drhov_vec[1]));
            G2 = -1.0 * u_roe * drho + (drhov_vec[0] * n[0] + drhov_vec[1] * n[1]);
            C1 = G1 * (s1 - fabs(lambda_3)) / pow(c_roe, 2) + G2 * s2 / c_roe;
            C2 = G1 * s2 / c_roe + (s1 - fabs(lambda_3)) * G2;

            ublas::matrix<double> FL = CalcAnalyticalFlux(uL, gamma);
            ublas::matrix<double> FR = CalcAnalyticalFlux(uR, gamma);

            ublas::vector<double> FL_hat(4, 0.0), FR_hat(4, 0.0);
            ublas::axpy_prod(FL, n, FL_hat, true);
            ublas::axpy_prod(FR, n, FR_hat, true);

            F_hat[0] = 0.5 * (FL_hat[0] + FR_hat[0]) - 0.5 * (fabs(lambda_3) * drho + C1);
            F_hat[1] = 0.5 * (FL_hat[1] + FR_hat[1]) - 0.5 * (fabs(lambda_3) * drhov_vec[0] + C1 * v_roe_vec[0] + C2 * n[0]);
            F_hat[2] = 0.5 * (FL_hat[2] + FR_hat[2]) - 0.5 * (fabs(lambda_3) * drhov_vec[1] + C1 * v_roe_vec[1] + C2 * n[1]);
            F_hat[3] = 0.5 * (FL_hat[3] + FR_hat[3]) - 0.5 * (fabs(lambda_3) * drhoE + C1 * H_roe + C2 * u_roe);

            mws = fabs(u_roe) + c_roe;

            return F_hat;

        }
        else
        {
            std::cout << "Unsupport flux name: " << type_flux << " Aborting" << std::endl;
            abort();
        }

    }

    ublas::vector<double> ApplyBoundaryCondition(ublas::vector<double> u, ublas::vector<double> norm,  char* boundary_type, Param& cparam, double &mws)
    {
        ublas::vector<double> num_flux(u.size());

        if (strcasecmp(boundary_type, "Inflow") == 0){
            double p = (cparam.gamma - 1.0) * (u[3] - 0.5 * (u[1]*u[1] + u[2]*u[2]) / u[0]);
            double un = (u[1] / u[0]) * norm[0] + (u[2] / u[0]) * norm[1];
            double c = sqrt(cparam.gamma * p / u[0]);
            double J = un + 2.0 * c / (cparam.gamma - 1); // Riemann Invariant
            double dn = cos(cparam.attack_angle)*norm[0] + sin(cparam.attack_angle)*norm[1];
            double R = 1.0;

            double Tt = 1.0 + 0.5 * (cparam.gamma - 1) * cparam.mach_inf * cparam.mach_inf;
            double pt = pow(Tt, cparam.gamma / (cparam.gamma - 1.0));

            // Solve for Mb
            // ca, cb, cc are coefficients for quadratic equation
            double tmpa = cparam.gamma * R * Tt * dn * dn - 0.5 * (cparam.gamma - 1.0) * J * J;
            double tmpb = 4.0 * cparam.gamma * R * Tt * dn / (cparam.gamma - 1.0);
            double tmpc = 4.0 * cparam.gamma * R * Tt  / pow((cparam.gamma - 1.0), 2) - J * J;
            double Mb1 = (-1.0 * tmpb - sqrt(tmpb * tmpb - 4.0 * tmpa * tmpc)) / (2 * tmpa);
            double Mb2 = (-1.0 * tmpb + sqrt(tmpb * tmpb - 4.0 * tmpa * tmpc)) / (2 * tmpa);
            double Mb;
            if (Mb1 < 0)
                Mb = Mb2;
            else
                Mb = Mb1;
            // Calculate the exterior states
            double Tb = Tt / (1.0 + 0.5*(cparam.gamma - 1.0) * Mb * Mb);
            double pb = pt * pow(Tb / Tt, cparam.gamma / (cparam.gamma - 1.0));
            double rhob = pb / (R * Tb);
            double cb = sqrt(cparam.gamma * pb / rhob);
            double vb[2] = {Mb * cb * cos(cparam.attack_angle), Mb * cb * sin(cparam.attack_angle)};
            double rhoEb = pb / (cparam.gamma - 1.0) + 0.5 * rhob * (vb[0] * vb[0] + vb[1] * vb[1]);
            ublas::vector<double> ub(4, 0.0);
            ub(0) = rhob; ub(1) = rhob * vb[0]; ub(2) = rhob * vb[1]; ub(3) = rhoEb;
            ublas::matrix<double> F_b = CalcAnalyticalFlux(ub, cparam.gamma);
            ublas::axpy_prod(F_b, norm, num_flux, true);
            // Max Wave Speed
            mws = fabs(vb[0] * norm[0] + vb[1] * norm[1]) + cb;

        } else if (strcasecmp(boundary_type, "Inviscid_Wall") == 0){

            double vb[2];

            vb[0] = u[1] / u[0] - ((u[1] / u[0]) * norm[0] + (u[2] / u[0]) * norm[1]) * norm[0];
            vb[1] = u[2] / u[0] - ((u[1] / u[0]) * norm[0] + (u[2] / u[0]) * norm[1]) * norm[1];

            double pb = (cparam.gamma - 1.0) * (u[3] - 0.5 * u[0] * (vb[0] * vb[0] + vb[1] * vb[1]));

            num_flux[0] = 0.0;
            num_flux[1] = norm[0] * pb;
            num_flux[2] = norm[1] * pb;
            num_flux[3] = 0.0;

            double p = (cparam.gamma - 1.0) * (u[3] - 0.5 * (u[1]*u[1] + u[2]*u[2]) / u[0]);
            mws = sqrt(pow(u[1] / u[0], 2) + pow(u[2] / u[0], 2)) + sqrt(cparam.gamma * p / u[0]);

        } else if (strcasecmp(boundary_type, "Subsonic_Outflow") == 0){

            /* Interior entropy*/
            double p = (cparam.gamma - 1.0) * (u[3] - 0.5 * (u[1]*u[1] + u[2]*u[2]) / u[0]);
            double S = p / pow(u[0], cparam.gamma);
            double pb = cparam.p_inf;
            double rhob = pow(pb / S,  1.0 / cparam.gamma);
            double cb = sqrt(cparam.gamma * pb / rhob);
            double un = (u[1] / u[0]) * norm[0] + (u[2] / u[0]) * norm[1];

            double c = sqrt(cparam.gamma * p / u[0]);
            double J = un + 2.0 * c / (cparam.gamma - 1.0); // Riemann Invariant

            double ub_n = J - 2.0 * cb / (cparam.gamma - 1.0);

            /* Solve for vb*/
            double vb[2] = {0.0, 0.0};
            vb[0] = u[1] / u[0] - norm[0] * (u[1]*norm[0] + u[2]*norm[1]) / u[0] + ub_n * norm[0];
            vb[1] = u[2] / u[0] - norm[1] * (u[1]*norm[0] + u[2]*norm[1]) / u[0] + ub_n * norm[1];

            double rhoEb = pb / (cparam.gamma - 1.0) + 0.5 * rhob * (vb[0] * vb[0] + vb[1] * vb[1]);

            ublas::vector<double> ub(4, 0.0);
            ub(0) = rhob; ub(1) = rhob * vb[0]; ub(2) = rhob * vb[1]; ub(3) = rhoEb;
            ublas::matrix<double> F_b = CalcAnalyticalFlux(ub, cparam.gamma);
            ublas::axpy_prod(F_b, norm, num_flux, true); // num_flux is passed out

            mws = sqrt(pow(ub[1] / ub[0], 2) + pow(ub[2] / ub[0], 2)) + sqrt(cparam.gamma * pb / ub[0]);

        } else if (strcasecmp(boundary_type, "Free_Stream") == 0){
            double mws_temp;
            ublas::vector<double> u_free(4);

            u_free[0] = 1.0;
            u_free[1] = cparam.mach_inf * cos(cparam.attack_angle);
            u_free[2] = cparam.mach_inf * sin(cparam.attack_angle);
            u_free[3] = 1 / (cparam.gamma * (cparam.gamma - 1)) +
                            0.5 * cparam.mach_inf * cparam.mach_inf;

            num_flux = CalcNumericalFlux(u, u_free, norm, cparam.gamma, "roe", mws_temp);
            mws = mws_temp;
        } else{
            printf("ERROR: Unknown Boundary Condition: %s \n", boundary_type);
            abort();
        }
        return num_flux;
    }
} // end namespace euler