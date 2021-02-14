#include <vector>
#include <math.h>

#include "Kinetics.hpp"
#include "TDMA_solver.hpp"
#include "Thermodynamic_Properties.hpp"

#define MAX_T_STEPS 100000
#define N 500

#define L 6.35E-3   // m
#define D 6.35E-3   // m

#define Dt 1E-5     // s
#define Dx 12.7E-6  // m

#define phi_P 0.5

#define T_a 298.0   // K
#define T_i 988.0   // K
#define T_f 1911.0  // K

vector<vector<float>> T(MAX_T_STEPS, vector<float>(N+1, T_a));

vector<float> E(N-1), F(N-1), G(N-1), B(N-1), T_Soln(N-1), eta(N+1, 0.0);
unsigned int R_Zone_Start_Index, R_Zone_Stop_Index, i;
float kDt;

TDMA_solver::solver my_solver(N-1);

float Q_dot_conv(float Temp);
float Q_dot_rad(float Temp);

int main(int argc, char const *argv[])
{
    float lambda_m_P    = get_lambda_m(lambda_p_P, lambda_Ar_P, phi_P, 0.5);
    float lambda_m_R    = get_lambda_m(lambda_p_R, lambda_Ar_R, phi_P, 0.5);
    float lambda_m_PC   = lambda_NiAl;

    float rho_m_P   = get_rho_m(phi_P, rho_Ar_P);
    float rho_m_R   = get_rho_m(phi_P, rho_Ar_R);
    float rho_m_PC  = rho_NiAl;

    float Cp_m_P    = get_Cp_m(phi_P, Cp_p_P, Cp_Ar_P, rho_Ar_P);
    float Cp_m_R    = get_Cp_m(phi_P, Cp_p_R, Cp_Ar_R, rho_Ar_R);
    float Cp_m_PC   = Cp_NiAl;
    
    T[0] = vector<float>(N+1, T_f);
    R_Zone_Start_Index  = 0;
    R_Zone_Stop_Index   = 1;
    
    for (int n=1; n<MAX_T_STEPS; n++)
    {
        // Post Combustion Zone
        for (i=0; i<R_Zone_Start_Index; i++)
        {
            E[i] = lambda_m_PC / pow(Dx,2);
            F[i] = - ( 2*lambda_m_PC/pow(Dx,2) + rho_m_PC*Cp_m_PC/Dt);
            G[i] = lambda_m_PC / pow(Dx,2);
            B[i] = Q_dot_conv(T[n-1][i+1]) + Q_dot_rad(T[n-1][i+1]) - rho_m_PC*Cp_m_PC*T[n-1][i+1]/Dt;
        }

        // Reaction Zone
        for (i=R_Zone_Start_Index; i<R_Zone_Stop_Index; i++)
        {
            E[i] = lambda_m_R / pow(Dx,2);
            F[i] = - ( 2*lambda_m_R/pow(Dx,2) + rho_m_R*Cp_m_R/Dt);
            G[i] = lambda_m_R / pow(Dx,2);
            B[i] = Q_dot_conv(T[n-1][i+1]) + Q_dot_rad(T[n-1][i+1]) -  - rho_m_R*Cp_m_R*T[n-1][i+1]/Dt;
        }

        // Preheat Zone
        for (i=R_Zone_Stop_Index; i<N-1; i++)
        {
            E[i] = lambda_m_P / pow(Dx,2);
            F[i] = - ( 2*lambda_m_P/pow(Dx,2) + rho_m_P*Cp_m_P/Dt);
            G[i] = lambda_m_P / pow(Dx,2);
            B[i] = Q_dot_conv(T[n-1][i+1]) + Q_dot_rad(T[n-1][i+1]) - rho_m_P*Cp_m_P*T[n-1][i+1]/Dt;
        }

        B[0]    -= E[0] * T_f;
        F[N-2]  += lambda_m_P/pow(Dx,2);
        
        // Setting up solver and solving
        my_solver.setup_banded_matrix(E, F, G);
        T_Soln = my_solver.solve(B);

        copy(T_Soln.begin(), T_Soln.end(), T[n].begin()+1);

        for (i=0; i<N; i++)
        {
            if (i-1 > R_Zone_Start_Index)
            {
                R_Zone_Start_Index = i-1;
            }

            kDt = Dt / get_t_b(T[n][i+1]);
        }

    }
    return 0;
}

float Q_dot_conv(float Temp)
{
    return 4 * 19.68 * (Temp - T_a) / D;
}

float Q_dot_rad(float Temp)
{
    return 4 * 0.25 * 5.6704E-8 * (pow(Temp, 4) - pow(Temp, 4)) / D;
}
