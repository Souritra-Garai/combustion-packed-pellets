#include "Temperature_Step_Iterator.hpp"

// Function definitions for class Temperature Iterator Base

Temperature_Iterator_Base::Temperature_Iterator_Base(unsigned int n)
{
    N = n;
    T_ign = 0;
}

void Temperature_Iterator_Base::reset_iterators()
{ 
    i = 0;
    reset_other_iterators();
}

bool Temperature_Iterator_Base::iterators_in_range() { return i < N; }

bool Temperature_Iterator_Base::in_reaction_zone(float T)
{
    return (T > T_ign && iterators_in_range());
}

bool Temperature_Iterator_Base::in_post_combustion_zone(float eta)
{
    return (eta > 0.99999 && iterators_in_range());
}

void Temperature_Iterator_Base::increment_iterators()
{ 
    i++;
    increment_other_iterators();
}

void Temperature_Iterator_Base::increment_iterators(vector<float>::iterator T)
{
    increment_iterators();
    T++;
}

void Temperature_Iterator_Base::increment_iterators(vector<float>::iterator T, vector<float>::iterator eta)
{
    increment_iterators(T);
    eta++;
}

void Temperature_Iterator_Base::increment_iterators(vector<float>::iterator T, vector<float>::iterator eta, vector<float>::iterator T_new)
{
    increment_iterators(T, eta);
    T_new++;
}

void Temperature_Iterator_Base::reset_other_iterators() {;}

void Temperature_Iterator_Base::increment_other_iterators() {;}

void Temperature_Iterator_Base::set_ignition_temperature(float T)
{
    T_ign = T;
}

// Function declarations for class Temperature Predictor Step Iterator

Temperature_Predictor_Iterator::Temperature_Predictor_Iterator(unsigned int n): 
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
    b_P = b_R = b_PC = 0;

    reset_iterators();
}

void Temperature_Predictor_Iterator::reset_other_iterators()
{
    E = E_arr.begin();
    F = F_arr.begin();
    G = G_arr.begin();
    B = B_arr.begin();
}

void Temperature_Predictor_Iterator::increment_other_iterators()
{
    E++;
    F++;
    G++;
    B++;
}

void Temperature_Predictor_Iterator::setup_banded_matrix(   
    vector<float>::iterator T, 
    vector<float>::iterator eta,
    void (*set_BC_X0)(float &e, float &f, float &g, float &b),
    void (*set_BC_XN)(float &e, float &f, float &g, float &b)    )
{
    for (reset_iterators(); in_post_combustion_zone(*eta); increment_iterators(T, eta))
    {
        *E = e_PC;
        *F = f_PC;
        *G = g_PC;
        *B = b_PC*(*T);
    }

    for (; in_reaction_zone(*T); increment_iterators(T))
    {
        *E = e_R;
        *F = f_R;
        *G = g_R;
        *B = b_R*(*T);
    }
    
    for (; iterators_in_range(); increment_iterators(T))
    {
        *E = e_P;
        *F = f_P;
        *G = g_P;
        *B = b_P*(*T);
    }

    set_BC_X0(E_arr[0], F_arr[0], G_arr[0], B_arr[0]);
    set_BC_XN(E_arr[N-1], F_arr[N-1], G_arr[N-1], B_arr[N-1]);

    my_solver.setup_banded_matrix(E_arr, F_arr, G_arr);
}

vector<float> Temperature_Predictor_Iterator::get_solution()
{
    return my_solver.solve(B_arr);
}

void Temperature_Predictor_Iterator::assign_coefficients_P  (float lambda, float rho, float Cp, float Dx, float Dt)
{
    e_P = calc_e(lambda, Dx);
    f_P = calc_f(lambda, rho, Cp, Dx, Dt);
    g_P = calc_g(lambda, Dx);
    b_P = calc_b(rho, Cp, Dt);
}

void Temperature_Predictor_Iterator::assign_coefficients_PC (float lambda, float rho, float Cp, float Dx, float Dt)
{
    e_PC = calc_e(lambda, Dx);
    f_PC = calc_f(lambda, rho, Cp, Dx, Dt);
    g_PC = calc_g(lambda, Dx);
    b_PC = calc_b(rho, Cp, Dt);
}

void Temperature_Predictor_Iterator::assign_coefficients_R  (float lambda, float rho, float Cp, float Dx, float Dt)
{
    e_R = calc_e(lambda, Dx);
    f_R = calc_f(lambda, rho, Cp, Dx, Dt);
    g_R = calc_g(lambda, Dx);
    b_R = calc_b(rho, Cp, Dt);
}

float Temperature_Predictor_Iterator::calc_e(float lambda, float Dx) {return lambda / (Dx * Dx);}
float Temperature_Predictor_Iterator::calc_g(float lambda, float Dx) {return lambda / (Dx * Dx);}
float Temperature_Predictor_Iterator::calc_b(float rho, float Cp, float Dt) {return rho * Cp / Dt;}
float Temperature_Predictor_Iterator::calc_f(float lambda, float rho, float Cp, float Dx, float Dt)
{
    return - (2 * lambda / (Dx * Dx) + rho * Cp / Dt);
}

// Function declarations for Temperature Corrector step iterator

Temperature_Corrector_Iterator::Temperature_Corrector_Iterator(unsigned int n) : Temperature_Iterator_Base(n)
{
    time_derivative_coeff_P = time_derivative_coeff_R = time_derivative_coeff_PC = 0;
    heat_source_term_coeff_R = 0;
    heat_conv_rad_const = 0;
    heat_rad_coeff_num = 0;
    heat_rad_coeff_den = 0;
    heat_conv_coeff = 0;
}

void Temperature_Corrector_Iterator::set_coeffs_pre_heat_zone(float rho, float Cp, float Delta_t)
{
    time_derivative_coeff_P = rho * Cp / Delta_t;
}

void Temperature_Corrector_Iterator::set_coeffs_reaction_zone(float rho, float Cp, float Delta_t)
{
    time_derivative_coeff_R = rho * Cp / Delta_t;
    heat_source_term_coeff_R = rho * Q_r / MW_Prod;
}

void Temperature_Corrector_Iterator::set_coeffs_post_combustion_zone(float rho, float Cp, float Delta_t)
{
    time_derivative_coeff_PC = rho * Cp / Delta_t;
}

void Temperature_Corrector_Iterator::set_conv_rad_properties(float h, float epsilon, float D, float T_a)
{
    heat_conv_coeff = 4 * h / D;
    heat_rad_coeff_num = 12 * epsilon * Stefan_Boltzmann_constant / D;
    heat_rad_coeff_den = 16 * epsilon * Stefan_Boltzmann_constant / D;
    heat_conv_rad_const = 4 * h * T_a / D + 4 * epsilon * Stefan_Boltzmann_constant * pow(T_a, 4) / D;
}

void Temperature_Corrector_Iterator::temperature_update(vector<float>::iterator T_hat, vector<float>::iterator eta, vector<float>::iterator T_new)
{
    reset_iterators();

    for (; in_post_combustion_zone(*eta); increment_iterators(T_hat, eta, T_new))
    {
        *T_new = (
            heat_rad_coeff_num * pow(*T_hat, 4) +
            time_derivative_coeff_PC * (*T_hat) +
            heat_conv_rad_const
        ) / (
            heat_rad_coeff_den * pow(*T_hat, 3) +
            time_derivative_coeff_PC +
            heat_conv_coeff
        );
    }

    for (; in_reaction_zone(*T_hat); increment_iterators(T_hat, eta, T_new))
    {
        *T_new = ( -
            heat_source_term_coeff_R * (1 - (*eta)) * (calc_k(*T_hat) - calc_k_prime(*T_hat)*(*T_hat)) +
            heat_rad_coeff_num * pow(*T_hat, 4) +
            time_derivative_coeff_R * (*T_hat) +
            heat_conv_rad_const
        ) / ( -
            heat_source_term_coeff_R * (1 - (*eta)) * calc_k_prime(*T_hat) +
            heat_rad_coeff_den * pow(*T_hat, 3) +
            time_derivative_coeff_R +
            heat_conv_coeff
        );
    }

    for (; iterators_in_range(); increment_iterators(T_hat, eta, T_new))
    {
        *T_new = (
            heat_rad_coeff_num * pow(*T_hat, 4) +
            time_derivative_coeff_P * (*T_hat) +
            heat_conv_rad_const
        ) / (
            heat_rad_coeff_den * pow(*T_hat, 3) +
            time_derivative_coeff_P +
            heat_conv_coeff
        );
    }
}
