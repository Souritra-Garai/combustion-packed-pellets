#include "Temperature_Iterator.hpp"

// Function declarations for class Temperature Predictor Step Iterator

Temperature_Iterator::Temperature_Iterator(unsigned int n, float del_t): 
    E_arr(n,0),
    F_arr(n,0),
    G_arr(n,0),
    B_arr(n,0), 
    my_solver(n),
    Temperature_Iterator_Base(n)
{
    Delta_t = del_t;

    e_P = e_R = e_PC = 0;
    f_0_P = f_0_R = f_0_PC = 0;
    g_P = g_R = g_PC = 0;

    b_1_P = b_1_R = b_1_PC = 0;
    b_0 = 0;

    reaction_term_coeff = 0;

    reset_iterators();
}

void Temperature_Iterator::reset_other_iterators()
{
    E = E_arr.begin();
    F = F_arr.begin();
    G = G_arr.begin();
    B = B_arr.begin();
}

void Temperature_Iterator::increment_other_iterators()
{
    E++;
    F++;
    G++;
    B++;
}

void Temperature_Iterator::setup_banded_matrix(   
    vector<float>::iterator T, 
    vector<float>::iterator eta,
    vector<float>::iterator omega,
    void (*set_BC_X0)(float &e, float &f, float &g, float &b),
    void (*set_BC_XN)(float &e, float &f, float &g, float &b)    )
{
    for (reset_iterators(); in_post_combustion_zone(*eta); increment_iterators(T, eta, omega))
    {
        *E = e_PC;
        *G = g_PC;

        *F = f_0_PC + f_3 * pow(*T, 3);

        *B = b_0 + b_1_PC * (*T) + b_4 * pow(*T, 4);
    }

    for (; in_reaction_zone(*T); increment_iterators(T, eta, omega))
    {
        *E = e_R;
        *G = g_R;

        pair<float, float> reaction_terms = reaction_term(*eta, *omega, *T, Delta_t);

        *F = f_0_R + f_3 * pow(*T, 3) - reaction_term_coeff * reaction_terms.second;

        *B = b_0 + b_1_R * (*T) + b_4 * pow(*T, 4) + reaction_term_coeff * reaction_terms.first;
    }
    
    for (; iterators_in_range(); increment_iterators(T))
    {
        *E = e_P;
        *G = g_P;

        *F = f_0_P + f_3 * pow(*T, 3);

        *B = b_0 + b_1_P * (*T) + b_4 * pow(*T, 4);
    }

    set_BC_X0(E_arr[0], F_arr[0], G_arr[0], B_arr[0]);
    set_BC_XN(E_arr[N-1], F_arr[N-1], G_arr[N-1], B_arr[N-1]);

    my_solver.setup_banded_matrix(E_arr, F_arr, G_arr);
}

vector<float> Temperature_Iterator::get_solution()
{
    return my_solver.solve(B_arr);
}

void Temperature_Iterator::assign_coefficients_P  (float lambda, float rho, float Cp, float Delta_x)
{
    e_P = calc_e(lambda, Delta_x);
    g_P = calc_g(lambda, Delta_x);

    f_0_P = calc_f(lambda, rho, Cp, Delta_x, Delta_t);

    b_1_P = calc_b_1(rho, Cp, Delta_t);
}

void Temperature_Iterator::assign_coefficients_PC (float lambda, float rho, float Cp, float Delta_x)
{
    e_PC = calc_e(lambda, Delta_x);
    g_PC = calc_g(lambda, Delta_x);

    f_0_PC = calc_f(lambda, rho, Cp, Delta_x, Delta_t);

    b_1_PC = calc_b_1(rho, Cp, Delta_t);
}

void Temperature_Iterator::assign_coefficients_R  (float lambda, float rho, float Cp, float Delta_x)
{
    e_R = calc_e(lambda, Delta_x);
    g_R = calc_g(lambda, Delta_x);

    f_0_R = calc_f(lambda, rho, Cp, Delta_x, Delta_t);

    b_1_R = calc_b_1(rho, Cp, Delta_t);
}

float Temperature_Iterator::calc_e(float lambda, float Delta_x) {return - lambda / (Delta_x * Delta_x);}
float Temperature_Iterator::calc_g(float lambda, float Delta_x) {return - lambda / (Delta_x * Delta_x);}

float Temperature_Iterator::calc_b_1(float rho, float Cp, float Delta_t) {return rho * Cp / Delta_t;}

float Temperature_Iterator::calc_f(float lambda, float rho, float Cp, float Delta_x, float Delta_t)
{
    return  (2 * lambda / (Delta_x * Delta_x)) + (rho * Cp / Delta_t);
}

void Temperature_Iterator::set_curved_surface_heat_loss(float Dia, float h_CSA, float epsi, float T_atmos)
{
    b_0 = 4 * (h_CSA * T_atmos + epsi * Stefan_Boltzmann_constant * pow(T_atmos, 4)) / Dia;

    b_4 = 4 * 3 * epsi * Stefan_Boltzmann_constant / Dia;

    f_3 = 4 * 4 * epsi * Stefan_Boltzmann_constant / Dia;

    f_0_P  += 4 * h_CSA / Dia;
    f_0_R  += 4 * h_CSA / Dia;
    f_0_PC += 4 * h_CSA / Dia;
}

void Temperature_Iterator::set_reaction_term_coeff(float density, float volume_fraction, float mol_wt)
{
    reaction_term_coeff = volume_fraction * density / mol_wt;
}