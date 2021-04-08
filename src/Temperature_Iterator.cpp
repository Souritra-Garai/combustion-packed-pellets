#include "Temperature_Iterator.hpp"

#include <iostream>
#include <iterator> // needed for std::ostram_iterator

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << '[';
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  return out;
}

Temperature_Iterator::Temperature_Iterator(unsigned int n) :
    T_VECTOR(n, 0), ETA_VECTOR(n, 0), OMEGA_VECTOR(n, 0),
    E_VECTOR(n, 0), F_VECTOR(n, 0), G_VECTOR(n, 0), B_VECTOR(n, 0),
    SOLVER(n)
{
    N = n;

    L = D = 0;

    T_atm = T_ign = T_final = 0;

    Delta_t = Delta_x = 0;

    rho = C_p = lambda = nullptr;

    molar_rho_limiting_reactant = 0;

    h_conv = epsilon = 0;
}

long double Temperature_Iterator::Calc_E(unsigned int i) {return - lambda[i] / (Delta_x * Delta_x);}
long double Temperature_Iterator::Calc_G(unsigned int i) {return - lambda[i] / (Delta_x * Delta_x);}
long double Temperature_Iterator::Calc_F(unsigned int i, long double Temperature)
{
    return (rho[i] * C_p[i] / Delta_t) + (2 * lambda[i] / (Delta_x * Delta_x)) + 16 * epsilon * Stefan_Boltzmann_Constant * pow(Temperature, 3) / D + 4 * h_conv / D;
}
long double Temperature_Iterator::Calc_B(unsigned int i, long double Temperature)
{
    return rho[i] * C_p[i] * Temperature / Delta_t + 4 * h_conv * T_atm / D + 4 * epsilon * Stefan_Boltzmann_Constant * (3 * pow(Temperature, 4) + pow(T_atm, 4)) / D;
}

void Temperature_Iterator::Apply_BC_Banded_Matrix()
{
    // Boundary Condition at x = 0
    F_VECTOR[0] = - lambda[2] / Delta_x - h_conv - 4 * epsilon * Stefan_Boltzmann_Constant * pow(T_VECTOR[0], 3);
    G_VECTOR[0] = lambda[2] / Delta_x;
    // F_VECTOR[0] = 1;
    // G_VECTOR[0] += E_VECTOR[0];

    // Boundary Condition at x = L
    F_VECTOR[N-1] += G_VECTOR[N-1] + 4 * epsilon * Stefan_Boltzmann_Constant * pow(T_VECTOR[N-1], 3) / Delta_x + h_conv / Delta_x;
}

void Temperature_Iterator::Apply_BC_B_Vector()
{
    // Boundary Condition at x = 0
    B_VECTOR[0] = - h_conv * T_atm - epsilon * Stefan_Boltzmann_Constant * (3 * pow(T_VECTOR[0], 4) + pow(T_atm, 4));
    // B_VECTOR[0] = T_final;

    // Boundary Condition at x = L
    B_VECTOR[N-1] += h_conv * T_atm / Delta_x + epsilon * Stefan_Boltzmann_Constant * (3 * pow(T_VECTOR[N-1], 4) + pow(T_atm, 4)) / Delta_x;
}

void Temperature_Iterator::Setup_Matrix_Equation()
{
    // std::cout << OMEGA_VECTOR << std::endl;
    // std::cout << ETA_VECTOR << std::endl;
    // std::cout << T_VECTOR << std::endl;

    int i = 0;

    for (   Reset_Banded_Matrix_Iterators(), B = B_VECTOR.begin(), T = T_VECTOR.begin(), ETA = ETA_VECTOR.begin(), OMEGA = OMEGA_VECTOR.begin();
            In_Range_Banded_Matrix_Iterators() && In_Post_Combustion_Zone(*ETA, *T);
            Increment_Banded_Matrix_Iterators(), B++, T++, ETA++, OMEGA++   )
    {
        *E = Calc_E(2);
        *G = Calc_G(2);
        *F = Calc_F(2, *T);
        *B = Calc_B(2, *T);
        i++;
    }

    std::cout << "\tReaction zone start grid #\t" << i << std::endl;

    for (   ;
            In_Range_Banded_Matrix_Iterators() && In_Reaction_Zone(*T);
            Increment_Banded_Matrix_Iterators(), B++, T++, ETA++, OMEGA++   )
    {
        *E = Calc_E(1);
        *G = Calc_G(1);

        std::pair<long double, long double> reaction_terms = Calc_Reaction_Terms(*ETA, *T, Delta_t);

        std::cout << "\t\tB\t" << Calc_B(1, *T) << "\tF\t" << Calc_F(1, *T) << std::endl;
        std::cout << "\tReaction Term\t" << reaction_terms.first * molar_rho_limiting_reactant << '\t' << reaction_terms.second * molar_rho_limiting_reactant << std::endl;
        std::cout << "\t\tCalc\t" << (1- (*ETA)) << std::endl;

        *F = Calc_F(1, *T) - molar_rho_limiting_reactant * reaction_terms.second;
        *B = Calc_B(1, *T) + molar_rho_limiting_reactant * reaction_terms.first;

        std::cout << "\t\tTemperature\t" << *T << "\tOmega\t" << *OMEGA << std::endl;
        i++;
    }

    std::cout << "\tReaction zone end grid #\t" << i << std::endl;

    for (   ;
            In_Range_Banded_Matrix_Iterators();
            Increment_Banded_Matrix_Iterators(), B++, T++, ETA++, OMEGA++   )
    {
        *E = Calc_E(0);
        *G = Calc_G(0);
        *F = Calc_F(0, *T);
        *B = Calc_B(0, *T);
    }

    Apply_BC_Banded_Matrix();
    Apply_BC_B_Vector();

    // std::cout << E_VECTOR << std::endl;
    // std::cout << F_VECTOR << std::endl;
    // std::cout << G_VECTOR << std::endl;
    // std::cout << B_VECTOR << std::endl;

    SOLVER.setup_banded_matrix(E_VECTOR, F_VECTOR, G_VECTOR);
}

bool Temperature_Iterator::In_Reaction_Zone(long double Temperature) {return Temperature >= T_ign;}

bool Temperature_Iterator::In_Post_Combustion_Zone(long double Conversion, long double Temperature) {return Temperature >= T_final;} // Conversion > 1 - 0.00001;}

void Temperature_Iterator::Set_Time_Step_Length(long double Dt) {Delta_t = Dt;}

void Temperature_Iterator::Set_Pellet_Dimensions(long double Length, long double Diameter)
{
    L = Length;
    D = Diameter;

    Delta_x = L / (N-1);
}

void Temperature_Iterator::Set_Thermophysical_Properties(long double* Density, long double* Heat_Capacity, long double* Heat_Diffusivity)
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

void Temperature_Iterator::Set_Temperatures(long double Atmospheric_Temperature, long double Final_Temperature, long double Ignition_Temperature)
{
    T_atm   = Atmospheric_Temperature;
    T_ign   = Ignition_Temperature;
    T_final = Final_Temperature;
}

void Temperature_Iterator::Apply_Initial_Condition(std::vector<long double> T_INITIAL_VECTOR)
{
    std::copy(T_INITIAL_VECTOR.begin(), T_INITIAL_VECTOR.end(), T_VECTOR.begin());

    ETA_VECTOR[0] = 1l;

    for (   T = T_VECTOR.begin(), ETA = ETA_VECTOR.begin(), OMEGA = OMEGA_VECTOR.begin();
            T < T_VECTOR.end();
            T++, ETA++, OMEGA++ )
    {
        *OMEGA = Calc_Omega(*ETA, *T);
    }
}

void Temperature_Iterator::Set_Curved_Surface_Heat_Losses(long double Convective_Heat_Transfer_Coefficient, long double Emissivity)
{
    h_conv = Convective_Heat_Transfer_Coefficient;
    epsilon = Emissivity;
}

void Temperature_Iterator::Update_Reaction_Zone()
{
    for (   T = T_VECTOR.begin(), ETA = ETA_VECTOR.begin(), OMEGA = OMEGA_VECTOR.begin();
            T < T_VECTOR.end() && In_Post_Combustion_Zone(*ETA, *T);
            T++, ETA++, OMEGA++ )   ;

    for (   ;
            T < T_VECTOR.end() && In_Reaction_Zone(*T);
            T++, ETA++, OMEGA++ )
    {
        *OMEGA  = Omega_Update  (*ETA, *T, Delta_t);
        *ETA    = Eta_Update    (*ETA, *OMEGA, Delta_t);
    }

}

void Temperature_Iterator::Set_Molar_Density_Limiting_Reactant(long double Density, long double Volume_Fraction, long double Molar_Weight)
{
    molar_rho_limiting_reactant = Density * Volume_Fraction / Molar_Weight;
}
