#ifndef __TEMP_PRED_ITER__
#define __TEMP_PRED_ITER__

#include <vector>

#include "TDMA_solver.hpp"

using namespace std;

class Temperature_Predictor_Iterator
{
    private:

        unsigned int N;

        float e_P, e_R, e_PC;
        float f_P, f_R, f_PC;
        float g_P, g_R, g_PC;
        float b_P, b_R, b_PC;

        float T_f, T_ign;

        vector<float> E_arr, F_arr, G_arr, B_arr;
        vector<float>::iterator E, F, G, B;

        TDMA_solver::solver my_solver;

        void reset_iterators();
        bool iterators_in_range();
        void increment_iterators();
        void increment_iterators(vector<float>::iterator);
        void increment_iterators(vector<float>::iterator, vector<float>::iterator);
        bool in_post_combustion_zone(float eta);
        bool in_reaction_zone(float T);

    public:

        Temperature_Predictor_Iterator  (unsigned int n);
        
        void assign_coefficients_P  (float e, float f, float g, float b);
        void assign_coefficients_R  (float e, float f, float g, float b);
        void assign_coefficients_PC (float e, float f, float g, float b);

        void set_temperatures(float T_Final, float T_Ignition);

        void Temperature_Predictor_Iterator::setup_banded_matrix(   
            vector<float>::iterator T, 
            vector<float>::iterator eta,
            void (*set_BC_X0)(float &e, float &f, float &g, float &b),
            void (*set_BC_XN)(float &e, float &f, float &g, float &b)    );

        vector<float> get_solution();
};


#endif