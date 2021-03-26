#include "Temperature_Iterator.hpp"

// Function declarations for class Temperature Predictor Step Iterator

Temperature_Iterator::Temperature_Iterator(unsigned int n): 
    E_arr(n,0),
    F_arr(n,0),
    G_arr(n,0),
    B_arr(n,0), 
    my_solver(n),
    Temperature_Iterator_Base(n)
{
    e_P = e_R = e_PC = 0;
    f_P = f_R = f_PC = 0;
    g_P = g_R = g_PC = 0;
    b_1_coef_P = b_1_coef_R = b_1_coef_PC = 0;
    b_static_P = b_static_R = b_static_PC = 0;
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
    for (reset_iterators(); in_post_combustion_zone(*eta); increment_iterators(T, eta))
    {
        *E = e_PC;
        *F = f_PC;
        *G = g_PC;
        *B = b_1_coef_PC*(*T);
    }

    for (; in_reaction_zone(*T); increment_iterators(T))
    {
        *E = e_R;
        *F = f_R;
        *G = g_R;
        *B = b_1_coef_R*(*T);
    }
    
    for (; iterators_in_range(); increment_iterators(T))
    {
        *E = e_P;
        *F = f_P;
        *G = g_P;
        *B = b_1_coef_P*(*T);
    }

    set_BC_X0(E_arr[0], F_arr[0], G_arr[0], B_arr[0]);
    set_BC_XN(E_arr[N-1], F_arr[N-1], G_arr[N-1], B_arr[N-1]);

    my_solver.setup_banded_matrix(E_arr, F_arr, G_arr);
}

vector<float> Temperature_Iterator::get_solution()
{
    return my_solver.solve(B_arr);
}

void Temperature_Iterator::assign_coefficients_P  (float lambda, float rho, float Cp, float Delta_x, float Delta_t)
{
    e_P = calc_e(lambda, Delta_x);
    f_P = calc_f(lambda, rho, Cp, D, h, Delta_x, Delta_t);
    g_P = calc_g(lambda, Delta_x);
    b_static_P = calc_b(rho, Cp, Delta_t);
}

void Temperature_Iterator::assign_coefficients_PC (float lambda, float rho, float Cp, float Delta_x, float Delta_t)
{
    e_PC = calc_e(lambda, Delta_x);
    f_PC = calc_f(lambda, rho, Cp, D, h, Delta_x, Delta_t);
    g_PC = calc_g(lambda, Delta_x);
    b_static_PC = calc_b(rho, Cp, Delta_t);
}

void Temperature_Iterator::assign_coefficients_R  (float lambda, float rho, float Cp, float Delta_x, float Delta_t)
{
    e_R = calc_e(lambda, Delta_x);
    f_R = calc_f(lambda, rho, Cp, D, h, Delta_x, Delta_t);
    g_R = calc_g(lambda, Delta_x);
    b_static_R = calc_b(rho, Cp, Delta_t);
}

float Temperature_Iterator::calc_e(float lambda, float Delta_x) {return - lambda / (Delta_x * Delta_x);}
float Temperature_Iterator::calc_g(float lambda, float Delta_x) {return - lambda / (Delta_x * Delta_x);}
float Temperature_Iterator::calc_b(float rho, float Cp, float Delta_t)
{
    return rho * Cp / Delta_t;
}
float Temperature_Iterator::calc_f(float lambda, float rho, float Cp, float D, float h, float Delta_x, float Delta_t)
{
    return 2 * lambda / (Delta_x * Delta_x) + rho * Cp / Delta_t - (4 / D)*h;
}
