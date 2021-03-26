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

        float Delta_t;

        float e_P, e_R, e_PC;
        float g_P, g_R, g_PC;
        float f_0_P, f_0_R, f_0_PC;
        float b_1_P, b_1_R, b_1_PC;
        
        float b_0, b_4, f_3;

        vector<float> E_arr, F_arr, G_arr, B_arr;
        vector<float>::iterator E, F, G, B;

        TDMA_solver::solver my_solver;

    protected:

        void reset_other_iterators();    
        void increment_other_iterators();
                
        float calc_e(float lambda, float Delta_x);
        float calc_g(float lambda, float Delta_x);
        float calc_f(float lambda, float rho, float Cp, float Delta_x, float Delta_t);
        float calc_b_1(float rho, float Cp, float Dt);
    
    public:

        Temperature_Iterator(unsigned int n, float delta_t);

        void assign_coefficients_P  (float lambda, float rho, float Cp, float Delta_x);
        void assign_coefficients_R  (float lambda, float rho, float Cp, float Delta_x);
        void assign_coefficients_PC (float lambda, float rho, float Cp, float Delta_x);

        void set_curved_surface_heat_loss(
            float diameter,
            float convective_heat_transfer_coeff,
            float emmisivity,
            float T_atmosphere  );

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