#ifndef __TEMP_ITER__
#define __TEMP_ITER__

#include <vector>
#include <math.h>

#include "Kinetics.hpp"
#include "TDMA_solver.hpp"
#include "Thermodynamic_Properties.hpp"
#include "Temperature_Iterator_Base.hpp"

using namespace std;

class Temperature_Iterator: public Temperature_Iterator_Base
{
    private:

        float D;
        float h;
        float epsilon;

        float e_P, e_R, e_PC;
        float f_P, f_R, f_PC;
        float g_P, g_R, g_PC;
        float b_1_coef_P, b_1_coef_R, b_1_coef_PC;
        float b_static_P, b_static_R, b_static_PC;

        vector<float> E_arr, F_arr, G_arr, B_arr;
        vector<float>::iterator E, F, G, B;

        TDMA_solver::solver my_solver;

    protected:

        void reset_other_iterators();    
        void increment_other_iterators();
                
        float calc_e(float lambda, float Delta_x);
        float calc_g(float lambda, float Delta_x);
        float calc_f(float lambda, float rho, float Cp, float D, float h, float Delta_x, float Delta_t);
        float calc_b(float rho, float Cp, float Dt);
    
    public:

        Temperature_Iterator(unsigned int n);

        void assign_coefficients_P  (float lambda, float rho, float Cp, float Delta_x, float Delta_t);
        void assign_coefficients_R  (float lambda, float rho, float Cp, float Delta_x, float Delta_t);
        void assign_coefficients_PC (float lambda, float rho, float Cp, float Delta_x, float Delta_t);

        void setup_banded_matrix(
            vector<float>::iterator T,
            vector<float>::iterator eta,
            vector<float>::iterator omega,
            void (*set_BC_X0)(float &e, float &f, float &g, float &b),
            void (*set_BC_XN)(float &e, float &f, float &g, float &b)
        );

        vector<float> get_solution();
};

#endif