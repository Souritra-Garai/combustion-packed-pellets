#include <vector>
#include <math.h>
#include <fstream>
#include <iostream>

#include "Kinetics.hpp"
#include "TDMA_solver.hpp"
#include "Thermodynamic_Properties.hpp"
#include "Temperature_Iterator.hpp"

using namespace std;

#define MAX_ITER_CONV   1000
#define MAX_T_STEPS     1000

#define N 100

#define L 6.35E-3   // m
#define D 6.35E-3   // m

#define Dt 1E-5     // s
#define Dx L/N

#define phi     0.5
#define alpha   0.5

#define T_a 298.0   // K
#define T_i 988.0   // K
#define T_f 1911.0  // K

#define FILENAME "solutions/T_Solution.csv"

// 2D Array to store Temperature evolution
// Axis 0 (index#1) is time point
// Axis 1 (index#2) is x point
vector<vector<float>> T_data(1, vector<float>(N+1, T_a));
// N grid points excluding the grid point at x=0
// and includeing the point at x=L

// Array for holding the conversion at a point
vector<float> eta_arr(N, 0.0);

// Everything for performing the temperature iteration

// Functions to set the boundary conditions
void set_BC_X0(float &e, float &f, float &g, float &b);
void set_BC_XN(float &e, float &f, float &g, float &b);

bool converged(vector<float>, vector<float>);

int main(int argc, char const *argv[])
{
    cout << "Starting program...\n";

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

    cout << "Thermophysical properties initialized...\n";

    T_data[0][0] = T_f;

    TPI.assign_coefficients_P  (lambda_m_P,     rho_m_P,    Cp_m_P,     Dx, Dt);
    TPI.assign_coefficients_PC (lambda_m_PC,    rho_m_PC,   Cp_m_PC,    Dx, Dt);
    TPI.assign_coefficients_R  (lambda_m_R,     rho_m_R,    Cp_m_R,     Dx, Dt);

    TPI.set_ignition_temperature(T_i);

    TCI.set_coeffs_post_combustion_zone(rho_m_PC, Cp_m_PC, Dt);
    TCI.set_coeffs_pre_heat_zone(rho_m_P, Cp_m_P, Dt);
    TCI.set_coeffs_reaction_zone(rho_m_R, Cp_m_R, Dt);
    
    TCI.set_conv_rad_properties(19.68 , 0.25, D, T_a);
    
    TCI.set_ignition_temperature(T_i);

    cout << "Temperature step iterators initialized...\n";

    for (int n=1; n < MAX_T_STEPS; n++)
    {        
        vector<float> T_nk, T_hat;
        T_nk = T_data[n-1];

        unsigned int i = 0;

        do
        {
            TPI.setup_banded_matrix(T_nk.begin()+1, eta_arr.begin(), &set_BC_X0, &set_BC_XN);
            T_hat = TPI.get_solution();

            update_eta(eta_arr, T_hat, Dt);

            TCI.temperature_update(T_hat.begin(), eta_arr.begin(), T_nk.begin()+1);

            cout << "\tPerformed " << ++i << " iterations under time step " << n << endl;

        } while (!converged(T_hat, T_nk) && i < MAX_ITER_CONV);

        T_data.push_back(T_nk);       

        cout << "Completed " << n << " time steps\n";

        if (converged(T_data[n-1], T_data[n])) break;
    }

    cout << "Saving to file...\n";

    ofstream file(FILENAME);

    for (auto T_row = T_data.begin(); T_row < T_data.end(); T_row++)
    {
        for (auto T = T_row->begin(); T < T_row->end(); T++)
        {
            file << to_string((*T)) << ',';
        }

        file << endl;        
    }

    file.close();
    cout << "Solution saved to file " << FILENAME << endl;

    return 0;
}

// Function definitions to set the boundary conditions
void set_BC_X0(float &e, float &f, float &g, float &b) { b -= e * T_f; }

void set_BC_XN(float &e, float &f, float &g, float &b) { e += g; }

bool converged(vector<float> T_prev, vector<float> T_new)
{
    bool flag = true;

    for (auto T_p = T_prev.begin(), T_n = T_new.begin(); T_p < T_prev.end() && flag; T_p++, T_n++)
    {
        if ((*T_p - *T_n > 0.001) || (*T_p - *T_n < -0.001)) flag = false;
    }
    return flag;
}