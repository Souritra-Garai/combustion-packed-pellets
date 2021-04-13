#include "Kinetics/Reaction.hpp"

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
    std::pair<long double, long double> partial_derivatives = Reaction_Kinetics.Get_Partial_Derivative_Reaction_Rate(Initial_Concentration * (1 - eta), T);

    long double r_A = Reaction_Kinetics.Get_Reaction_Rate(Initial_Concentration * (1 - eta), T);

    long double coeff = - Enthalpy_of_Reaction * Particle_Volume_Fraction;

    return std::pair<long double, long double> (
        coeff * partial_derivatives.second / (1 - partial_derivatives.first * Dt),
        - coeff * (r_A - T * partial_derivatives.second) / (1 - partial_derivatives.first * Dt)
    );
}

void Reaction::Update_Conversion(long double &eta, long double T, long double Dt)
{
    eta += - Reaction_Kinetics.Get_Reaction_Rate(Initial_Concentration * (1 - eta), T) * Dt / Initial_Concentration;
}