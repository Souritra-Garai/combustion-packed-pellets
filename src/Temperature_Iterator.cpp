#include "Temperature_Iterator.hpp"

#include <iostream>

Temperature_Iterator::Temperature_Iterator(unsigned int n) :
    T_VECTOR(n, 0),
    E_VECTOR(n, 0), F_VECTOR(n, 0), G_VECTOR(n, 0), B_VECTOR(n, 0),
    SOLVER(n)
{
    N = n;

    L = D = 0;

    T_atm = T_final = 0;

    Delta_t = Delta_x = 0;

    rho = C_p = lambda = 0;

    h_conv = epsilon = 0;
    
    E = E_VECTOR.begin();
    F = F_VECTOR.begin();
    G = G_VECTOR.begin();
    B = B_VECTOR.begin();

    T = T_VECTOR.begin();
}

long double Temperature_Iterator::Calc_E() {return - lambda / (Delta_x * Delta_x);}
long double Temperature_Iterator::Calc_G() {return - lambda / (Delta_x * Delta_x);}
long double Temperature_Iterator::Calc_F(long double Temperature)
{
    return (rho / Delta_t) + (2 * lambda / (Delta_x * Delta_x)) + 16 * epsilon * Stefan_Boltzmann_Constant * pow(Temperature, 3) / D + 4 * h_conv / D;
}
long double Temperature_Iterator::Calc_B(long double Temperature)
{
    return rho * Temperature / Delta_t + 4 * h_conv * T_atm / D + 4 * epsilon * Stefan_Boltzmann_Constant * (3 * pow(Temperature, 4) + pow(T_atm, 4)) / D;
}

void Temperature_Iterator::Apply_BC_Banded_Matrix()
{
    // Boundary Condition at x = 0
    F_VECTOR[0] += E_VECTOR[0] + 4 * epsilon * Stefan_Boltzmann_Constant * pow(T_VECTOR[0], 3) / Delta_x + h_conv / Delta_x;

    // Boundary Condition at x = L
    F_VECTOR[N-1] += G_VECTOR[N-1] + 4 * epsilon * Stefan_Boltzmann_Constant * pow(T_VECTOR[N-1], 3) / Delta_x + h_conv / Delta_x;
}

void Temperature_Iterator::Apply_BC_B_Vector()
{
    // Boundary Condition at x = 0
    B_VECTOR[0] += h_conv * T_atm / Delta_x + epsilon * Stefan_Boltzmann_Constant * (3 * pow(T_VECTOR[0], 4) + pow(T_atm, 4)) / Delta_x;

    // Boundary Condition at x = L
    B_VECTOR[N-1] += h_conv * T_atm / Delta_x + epsilon * Stefan_Boltzmann_Constant * (3 * pow(T_VECTOR[N-1], 4) + pow(T_atm, 4)) / Delta_x;
}

void Temperature_Iterator::Setup_Matrix_Equation()
{
    for (   Reset_Banded_Matrix_Iterators(), B = B_VECTOR.begin(), T = T_VECTOR.begin();
            In_Range_Banded_Matrix_Iterators();
            Increment_Banded_Matrix_Iterators(), B++, T++   )
    {
        *E = Calc_E();
        *G = Calc_G();
        *F = Calc_F(*T);
        *B = Calc_B(*T);
    }

    Apply_BC_Banded_Matrix();
    Apply_BC_B_Vector();

    SOLVER.setup_banded_matrix(E_VECTOR, F_VECTOR, G_VECTOR);
}

void Temperature_Iterator::Set_Time_Step_Length(long double Dt) {Delta_t = Dt;}

void Temperature_Iterator::Set_Pellet_Dimensions(long double Length, long double Diameter)
{
    L = Length;
    D = Diameter;

    Delta_x = L / (N-1);
}

void Temperature_Iterator::Set_Thermophysical_Properties(long double Density, long double Heat_Capacity, long double Heat_Diffusivity)
{
    rho     = Density;
    C_p     = Heat_Capacity;
    lambda  = Heat_Diffusivity;
}

void Temperature_Iterator::Reset_Banded_Matrix_Iterators()
{
    E = E_VECTOR.begin();
    F = F_VECTOR.begin();
    G = G_VECTOR.begin();
}

bool Temperature_Iterator::In_Range_Banded_Matrix_Iterators() {return E < E_VECTOR.end();}

void Temperature_Iterator::Increment_Banded_Matrix_Iterators() {E++; F++; G++;}

std::vector<long double> Temperature_Iterator::Get_Solution()
{
    std::vector<long double> New_T_Vector = SOLVER.solve(B_VECTOR);

    std::copy(New_T_Vector.begin(), New_T_Vector.end(), T_VECTOR.begin());

    return New_T_Vector;
}

void Temperature_Iterator::Set_Temperatures(long double Atmospheric_Temperature, long double Final_Temperature)
{
    T_atm   = Atmospheric_Temperature;
    T_final = Final_Temperature;
}

void Temperature_Iterator::Apply_Initial_Condition(std::vector<long double> T_INITIAL_VECTOR)
{
    std::copy(T_INITIAL_VECTOR.begin(), T_INITIAL_VECTOR.end(), T_VECTOR.begin());
}

void Temperature_Iterator::Set_Curved_Surface_Heat_Losses(long double Convective_Heat_Transfer_Coefficient, long double Emissivity)
{
    h_conv = Convective_Heat_Transfer_Coefficient;
    epsilon = Emissivity;
}
