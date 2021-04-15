#include "Kinetics/Kinetics.hpp"

#include <iostream>

Kinetics::Kinetics(
    std::string N,
    long double A,
    long double E_a,
    long double n 
) : Name(N)
{
    Pre_Exponential_Factor = A;
    Activation_Energy = E_a;
    Order = n;
}

long double Kinetics::Get_Reaction_Rate(
    long double C,
    long double T
)
{
    return - Get_Rate_Constant(T) * pow(C, Order);
}

long double Kinetics::Get_Rate_Constant(long double T)
{
    return Pre_Exponential_Factor * exp(- Activation_Energy / (Universal_Gas_Constant * T));
}

std::pair<long double, long double> Kinetics::Get_Partial_Derivative_Reaction_Rate(
    long double C,
    long double T
)
{
    long double k = Get_Rate_Constant(T);

    return std::pair<long double, long double> (
        - Order * k * pow(C, Order - 1),
        - k * pow(C, Order) * Activation_Energy / (Universal_Gas_Constant * T * T)
    );
}
