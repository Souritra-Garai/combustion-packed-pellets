#include "Temperature_Predictor_Iterator.hpp"

Temperature_Predictor_Iterator::Temperature_Predictor_Iterator(unsigned int n): 
    N(n),
    E_arr(n,0),
    F_arr(n,0),
    G_arr(n,0),
    B_arr(n,0), 
    my_solver(n)
{
    e_P = e_R = e_PC = 0;
    f_P = f_R = f_PC = 0;
    g_P = g_R = g_PC = 0;
    b_P = b_R = b_PC = 0;

    reset_iterators();
}

void Temperature_Predictor_Iterator::reset_iterators()
{
    E = E_arr.begin();
    F = F_arr.begin();
    G = G_arr.begin();
    B = B_arr.begin();
}

void Temperature_Predictor_Iterator::increment_iterators()
{
    E++;
    F++;
    G++;
    B++;
}

void Temperature_Predictor_Iterator::increment_iterators(vector<float>::iterator T)
{
    increment_iterators();
    T++;
}

void Temperature_Predictor_Iterator::increment_iterators(vector<float>::iterator T, vector<float>::iterator eta)
{
    increment_iterators(T);
    eta++;
}

bool Temperature_Predictor_Iterator::iterators_in_range()
{
    return B < B_arr.end();
}

bool Temperature_Predictor_Iterator::in_post_combustion_zone(float eta)
{
    return (eta > 0.99999 && iterators_in_range());
}

bool Temperature_Predictor_Iterator::in_reaction_zone(float T)
{
    return (T > T_ign && iterators_in_range());
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

