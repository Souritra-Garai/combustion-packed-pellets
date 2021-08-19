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
        long double Order_A;
        long double Order_B;

    public:

        Kinetics(
            std::string Name,
            long double Pre_Exponential_Factor,
            long double Activation_Energy,
            long double Order_A,
            long double Order_B
        );

        long double Get_Rate_Constant(long double Temperature);
        long double Get_Conversion_Rate(long double Conversion, long double Temperature);
        std::pair<long double, long double> Get_Partial_Derivative_Conversion_Rate(long double Conversion, long double Temperature);
};

#endif