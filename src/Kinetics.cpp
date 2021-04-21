#include "Kinetics/Kinetics.hpp"

#include <iostream>

Kinetics::Kinetics(
    std::string N,
    long double A,
    long double E_a,
    long double m,
    long double n
) : Name(N)
{
    Pre_Exponential_Factor = A;
    Activation_Energy = E_a;
    Order_A = m;
    Order_B = n;
}

long double Kinetics::Get_Conversion_Rate(
    long double eta,
    long double T
)
{
    return Get_Rate_Constant(T) * pow(1 - eta, Order_A) * pow(eta, Order_B);
}

long double Kinetics::Get_Rate_Constant(long double T)
{
    return Pre_Exponential_Factor * exp(- Activation_Energy / (Universal_Gas_Constant * T));
}

std::pair<long double, long double> Kinetics::Get_Partial_Derivative_Conversion_Rate(
    long double eta,
    long double T
)
{
    long double k = Get_Rate_Constant(T);

    return std::pair<long double, long double> (
        - Order_A * k * pow(1 - eta, Order_A - 1) * pow(eta, Order_B) + Order_B * k * pow(1 - eta, Order_A) * pow(eta, Order_B - 1),
        k * pow(1 - eta, Order_A) * pow(eta, Order_B) * Activation_Energy / (Universal_Gas_Constant * T * T)
    );
}
