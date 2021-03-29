#include <iostream>
#include <fstream>
#include <vector>

#include "Temperature_Iterator.hpp"
#include "Thermodynamic_Properties.hpp"

#define MAX_TIME_ITER   10000

#define N   1000

#define L   1.0 // m
#define D   1.0 // m

#define Dt  1E-5    // s
#define Dx  L/N     // m

#define T_atm   298.0   // K
#define T_f     1911.0  // K

#define phi     0.5
#define alpha   0.5

// Matrix for Temperature at the grid points
// for each time point
std::vector<std::vector<long double>> T_MATRIX(1, std::vector<long double> (N, T_atm));

Temperature_Iterator TI(N);

bool has_converged(std::vector<long double>, std::vector<long double>);

int main(int argc, char** argv)
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
    
    std::cout << "Setting up solver..." << std::endl;

    T_MATRIX[0][0] = T_f;

    TI.Set_Temperatures(T_atm, T_f);
    TI.Set_Time_Step_Length(Dt);
    TI.Set_Pellet_Dimensions(L, D);
    TI.Set_Thermophysical_Properties(rho_m_R, Cp_m_R, lambda_m_R);
    TI.Set_Curved_Surface_Heat_Losses(19.68, 0.25);

    TI.Apply_Initial_Condition(T_MATRIX[0]);

    std::cout << "Solver set up." << std::endl;

    int n = 0;

    do
    {
        std::cout << "Time step : " << (n+1) * Dt << " s" << std::endl;

        TI.Setup_Matrix_Equation();
        
        auto T_VECTOR = TI.Get_Solution();

        T_MATRIX.push_back(T_VECTOR);

        n++;

        std::cout << "Completed " << n << " iterations\n" << std::endl;

    } while (n < MAX_TIME_ITER && !has_converged(T_MATRIX[n], T_MATRIX[n-1]));

    std::cout << "Saving to file...\n";

    std::ofstream file("solutions/solution_1.csv");

    for (auto T_VECTOR = T_MATRIX.begin(); T_VECTOR < T_MATRIX.end(); T_VECTOR++)
    {
        for (auto T = T_VECTOR->begin(); T < T_VECTOR->end(); T++)
        {
            file << std::to_string((*T)) << ',';
        }

        file << std::endl;
    }

    file.close();
    std::cout << "Solution saved to file.." << std::endl;

    return 0;    
}

bool is_close(long double a, long double b)
{
    if (-0.01 < (a - b) && (a - b) < 0.01) return true;

    else return false;
}

bool has_converged(std::vector<long double> VECTOR_A, std::vector<long double> VECTOR_B)
{
    bool flag = true;
    std::vector<long double>::iterator A, B;
    
    for (A = VECTOR_A.begin(), B = VECTOR_B.begin(); A < VECTOR_A.end() && flag; A++, B++)

        flag = is_close(*A, *B);

    return flag;
}