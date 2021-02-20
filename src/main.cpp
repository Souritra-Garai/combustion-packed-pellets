#include <vector>
#include <math.h>

#include "Kinetics.hpp"
#include "TDMA_solver.hpp"
#include "Thermodynamic_Properties.hpp"

#define MAX_ITER_CONV   100
#define MAX_T_STEPS     1000

#define N 500

#define L 6.35E-3   // m
#define D 6.35E-3   // m

#define Dt 1E-3     // s
#define Dx L/N

#define phi     0.5
#define alpha   0.5

#define T_a 298.0   // K
#define T_i 988.0   // K
#define T_f 1911.0  // K

// 2D Array to store Temperature evolution
// Axis 0 (index#1) is time point
// Axis 1 (index#2) is x point
vector<vector<float>> T_data(MAX_T_STEPS, vector<float>(N+1, T_a));
// N+2 grid points including those at boundary

// Arrays for holding coefficients of discretized PDE
// Ref. Section 2.2.1
vector<float> E(N), F(N), G(N);
// Arrays for holding RHS values of discretized PDE
vector<float> B(N);
// Solver for discretized PDEs 
TDMA_solver::solver my_solver(N);

// Array for holding the conversion at a point
vector<float> eta_arr(N, 0.0);

vector<float>::iterator T, e, f, g, b, eta;

bool is_close(float, float);

int main(int argc, char const *argv[])
{
    // Finding the mean thermophysical properties for diffrent zones
    // Pre heat zone
    float lambda_m_P    = calc_lambda_m(phi, alpha, lambda_p_P, lambda_Ar_P);
    float lambda_m_PC   = calc_lambda_m(phi, alpha, lambda_NiAl_PC, lambda_Ar_PC);
    float lambda_m_R    = calc_lambda_m(phi, alpha, 0.5*(lambda_p_R + lambda_NiAl_R), lambda_Ar_R);
    // Post combustion zone
    float rho_m_P   = calc_rho_m(phi, rho_p_P, rho_Ar_P);
    float rho_m_PC  = calc_rho_m(phi, rho_NiAl_PC, rho_Ar_PC);
    float rho_m_R   = calc_rho_m(phi, 0.5*(rho_p_R + rho_NiAl_R), rho_Ar_R);
    // Reaction zone
    float Cp_m_P    = calc_Cp_m(phi, rho_p_P, rho_Ar_P, Cp_p_P, Cp_Ar_P);
    float Cp_m_R    = calc_Cp_m(phi, 0.5*(rho_p_R + rho_NiAl_R), rho_Ar_R, 0.5*(Cp_p_R + Cp_NiAl_R), Cp_Ar_R);
    float Cp_m_PC   = calc_Cp_m(phi, rho_NiAl_PC, rho_Ar_PC, Cp_NiAl_PC, Cp_Ar_PC);
    
    // Setting the initial conditions on temperature
    T_data[0][0] = T_f;
    fill(T_data[0].begin()+1, T_data[0].end(), T_a);

    T = T_data[0].begin()+1;
    e = E.begin();
    f = F.begin();
    g = G.begin();
    b = B.begin();
    eta = eta_arr.begin();

    for (; is_close(*eta, 1) && eta < eta_arr.end(); T++, e++, f++, g++, b++, eta++)
    {
        *e = lambda_m_PC / (Dx*Dx);
        *f = - ( 2 * lambda_m_PC / (Dx*Dx) + rho_m_PC * Cp_m_PC / Dt );
        *g = lambda_m_PC / (Dx*Dx);
        *b = - rho_m_PC * Cp_m_PC * (*T) / Dt;
    }

    for (; *T > T_i && eta < eta_arr.end(); T++, e++, f++, g++, b++, eta++)
    {
        *e = lambda_m_R / (Dx*Dx);
        *f = - ( 2 * lambda_m_R / (Dx*Dx) + rho_m_R * Cp_m_R / Dt );
        *g = lambda_m_R / (Dx*Dx);
        *b = - rho_m_R * Cp_m_R * (*T) / Dt;
    }

    for (; eta < eta_arr.end(); T++, e++, f++, g++, b++, eta++)
    {
        *e = lambda_m_P / (Dx*Dx);
        *f = - ( 2 * lambda_m_P / (Dx*Dx) + rho_m_P * Cp_m_P / Dt );
        *g = lambda_m_P / (Dx*Dx);
        *b = - rho_m_P * Cp_m_P * (*T) / Dt;
    }
    
    for (int n=1; n<MAX_T_STEPS; n++)
    {
        vector<float> T_hat(N);

        my_solver.setup_banded_matrix(E, F, G);
        T_hat = my_solver.solve(B);

        update_eta(eta_arr, T_hat, Dt);
    }

    return 0;
}