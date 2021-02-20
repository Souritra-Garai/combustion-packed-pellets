#ifndef __TEMP_PRED_ITER__
#define __TEMP_PRED_ITER__

#include <vector>
#include <math.h>

#include "TDMA_solver.hpp"
#include "Kinetics.hpp"
#include "Thermodynamic_Properties.hpp"

using namespace std;

class Temperature_Iterator_Base
{
    protected:

        unsigned int N;
        
        float T_ign;
        unsigned int i;

        void reset_iterators();
        
        bool iterators_in_range();
        bool in_reaction_zone(float T);
        bool in_post_combustion_zone(float eta);

        void increment_iterators();
        void increment_iterators(vector<float>::iterator);
        void increment_iterators(vector<float>::iterator, vector<float>::iterator);
        void increment_iterators(vector<float>::iterator, vector<float>::iterator, vector<float>::iterator);

        virtual void reset_other_iterators();
        virtual void increment_other_iterators();

    public:

        Temperature_Iterator_Base(unsigned int n);
        void set_ignition_temperature(float T_Ignition);
};

class Temperature_Predictor_Iterator: public Temperature_Iterator_Base
{
    private:

        float e_P, e_R, e_PC;
        float f_P, f_R, f_PC;
        float g_P, g_R, g_PC;
        float b_P, b_R, b_PC;

        vector<float> E_arr, F_arr, G_arr, B_arr;
        vector<float>::iterator E, F, G, B;

        TDMA_solver::solver my_solver;

    protected:

        void reset_other_iterators();    
        void increment_other_iterators();
                
        float calc_e(float lambda, float Delta_x);
        float calc_f(float lambda, float rho, float Cp, float Delta_x, float Delta_t);
        float calc_g(float lambda, float Delta_x);
        float calc_b(float rho, float Cp, float Dt);

    public:

        Temperature_Predictor_Iterator(unsigned int n);
        
        void assign_coefficients_P  (float lambda, float rho, float Cp, float Delta_x, float Delta_t);
        void assign_coefficients_R  (float lambda, float rho, float Cp, float Delta_x, float Delta_t);
        void assign_coefficients_PC (float lambda, float rho, float Cp, float Delta_x, float Delta_t);

        void setup_banded_matrix(   
            vector<float>::iterator T, 
            vector<float>::iterator eta,
            void (*set_BC_X0)(float &e, float &f, float &g, float &b),
            void (*set_BC_XN)(float &e, float &f, float &g, float &b)    );

        vector<float> get_solution();
};

class Temperature_Corrector_Iterator: public Temperature_Iterator_Base
{
    private:

        float time_derivative_coeff_P, time_derivative_coeff_R, time_derivative_coeff_PC;
        float heat_source_term_coeff_R;
        float heat_conv_rad_const;
        float heat_rad_coeff_num;
        float heat_rad_coeff_den;
        float heat_conv_coeff;
        
    public:

        Temperature_Corrector_Iterator(unsigned int n);
        
        void set_coeffs_pre_heat_zone(float rho, float Cp, float Delta_t);
        void set_coeffs_reaction_zone(float rho, float Cp, float Delta_t);
        void set_coeffs_post_combustion_zone(float rho, float Cp, float Delta_t);

        void set_conv_rad_properties(float h, float epsilon, float D, float T_a);
        void temperature_update(vector<float>::iterator T_hat, vector<float>::iterator eta, vector<float>::iterator T_new);
};

#endif