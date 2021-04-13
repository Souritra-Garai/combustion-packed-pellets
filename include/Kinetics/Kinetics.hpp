#ifndef __KINETICS__
#define __KINETICS__

#include <math.h>

#include <string>
#include <vector>

#define Universal_Gas_Constant 8.314

class Kinetics
{
    private:

        std::string Name;

        long double Pre_Exponential_Factor;
        long double Activation_Energy;
        long double Order;    

    public:

        Kinetics(
            std::string Name,
            long double Pre_Exponential_Factor,
            long double Activation_Energy,
            long double Order
        );

        long double Get_Rate_Constant(long double Temperature);
        long double Get_Reaction_Rate(long double Concentration, long double Temperature);
        std::pair<long double, long double> Get_Partial_Derivative_Reaction_Rate(long double Concentration, long double Temperature);
};

#endif