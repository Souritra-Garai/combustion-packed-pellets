#include "Kinetics/Reaction.hpp"

// #include <iostream>

Reaction::Reaction(Kinetics K, long double Delta_H_r, long double C_A0, long double phi_p) : Reaction_Kinetics(K)
{
    Enthalpy_of_Reaction = Delta_H_r;
    Initial_Concentration = C_A0;
    Particle_Volume_Fraction = phi_p;
}

std::pair<long double, long double> Reaction::Get_Linear_Expression(
    long double T,
    long double eta,
    long double Dt
)
{
    std::pair<long double, long double> partial_derivatives = Reaction_Kinetics.Get_Partial_Derivative_Conversion_Rate(eta, T);

    long double omega_prev = Reaction_Kinetics.Get_Conversion_Rate(eta, T);

    long double mDH_phi_CA0 = - Enthalpy_of_Reaction * Particle_Volume_Fraction * Initial_Concentration;

    long double gamma = (omega_prev - T * partial_derivatives.second) / (1 - Dt * partial_derivatives.first);

    long double kappa = partial_derivatives.second / (1 - Dt * partial_derivatives.first);

    // std::cout << "Partial Derivatives\t" << partial_derivatives.first << '\t' << partial_derivatives.second << std::endl;

    // std::cout << "Initial Concentration\t" << Initial_Concentration << std::endl;

    return std::pair<long double, long double> (
        mDH_phi_CA0 * gamma,
        mDH_phi_CA0 * kappa
    );
}

void Reaction::Update_Conversion(long double &eta, long double T, long double Dt)
{
    eta += Reaction_Kinetics.Get_Conversion_Rate(eta, T) * Dt;
}