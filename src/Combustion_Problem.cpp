#include "Combustion_Problem.hpp"
#include <iostream>

// #include <iterator> // needed for std::ostram_iterator

// template <typename T>
// std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
//   if ( !v.empty() ) {
//     out << '[';
//     std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
//     out << "\b\b]";
//   }
//   return out;
// }

Combustion_Problem::Combustion_Problem(
    unsigned int n,
    long double Dt,
    Combustion_Pellet P,
    Reaction R,
    long double T_ign,
    long double T_ad
) : T_VECTOR(n, 0), ETA_VECTOR(n, 0),
    E_VECTOR(n, 0), F_VECTOR(n, 0), G_VECTOR(n, 0), B_VECTOR(n, 0),
    SOLVER(n), Pellet(P), Combustion_Reaction(R)
{
    N = n;
    
    Delta_t = Dt;

    Ignition_Temperature = T_ign;
    Adiabatic_Combustion_Temperature = T_ad;

    Delta_x = Pellet.Get_Pellet_Length() / (N - 1);
}

void Combustion_Problem::Reset_Equation_Iterators()
{
    E = E_VECTOR.begin();
    F = F_VECTOR.begin();
    G = G_VECTOR.begin();
    B = B_VECTOR.begin();

    T = T_VECTOR.begin();
    ETA = ETA_VECTOR.begin();
}

bool Combustion_Problem::In_Range_Equation_Iterators()
{
    return E < (E_VECTOR.end()-1);
}

void Combustion_Problem::Increment_Equation_Iterators()
{
    E++;
    F++;
    G++;
    B++;
    
    T++;
    ETA++;
}

bool Combustion_Problem::In_Post_Combustion_Zone(long double Conversion, long double Temperature)
{
    return (Conversion > 1.0l - 0.00001l) && In_Range_Equation_Iterators();
    // return (Temperature >= Adiabatic_Combustion_Temperature) && In_Range_Equation_Iterators();
}

bool Combustion_Problem::In_Reaction_Zone(long double Temperature)
{
    return (Temperature >= Ignition_Temperature) && In_Range_Equation_Iterators();
}

void Combustion_Problem::Setup_Solver()
{
    Reset_Equation_Iterators();

    Pellet.Setup_X0_Ambient_Heat_Loss_BC_Equation(*E, *F, *G, *B, *T, Delta_x, Delta_t);

    Increment_Equation_Iterators();

    for ( ; In_Post_Combustion_Zone(*ETA, *T); Increment_Equation_Iterators())

        Pellet.Setup_Post_Combustion_Zone_Equation(*E, *F, *G, *B, *T, Delta_x, Delta_t);

    for ( ; In_Reaction_Zone(*T); Increment_Equation_Iterators())
    {
        Pellet.Setup_Reaction_Zone_Equation(*E, *F, *G, *B, *T, Delta_x, Delta_t);

        std::pair<long double, long double> Reaction_Terms = Combustion_Reaction.Get_Linear_Expression(*T, *ETA, Delta_t);

        // std::cout << "In Reaction Zone" << std::endl;
        // std::cout << Reaction_Terms.first << '\t' << Reaction_Terms.second << std::endl;
        // std::cout << *F << '\t' << *B << std::endl;

        *F += - Reaction_Terms.second;
        *B += Reaction_Terms.first;
    }

    for ( ; In_Range_Equation_Iterators(); Increment_Equation_Iterators())
    
        Pellet.Setup_Pre_Heat_Zone_Equation(*E, *F, *G, *B, *T, Delta_x, Delta_t);

    int zone = In_Reaction_Zone(*T) ? (In_Post_Combustion_Zone(*ETA, *T) ? 2 : 1) : 0;
    Pellet.Setup_XM_Ambient_Heat_Loss_BC_Equation(*E, *F, *G, *B, *T, Delta_x, Delta_t, zone);

    // std::cout << E_VECTOR << std::endl;
    // std::cout << F_VECTOR << std::endl;
    // std::cout << G_VECTOR << std::endl;
    // std::cout << B_VECTOR << std::endl;

    SOLVER.setup_banded_matrix(E_VECTOR, F_VECTOR, G_VECTOR);
}

void Combustion_Problem::Solve_and_Update_State()
{
    T_VECTOR = SOLVER.solve(B_VECTOR);

    // std::cout << T_VECTOR << std::endl;

    Reset_Equation_Iterators();

    for ( ; In_Post_Combustion_Zone(*ETA, *T); Increment_Equation_Iterators()) ;

    for ( ; In_Reaction_Zone(*T); Increment_Equation_Iterators())
    {
        Combustion_Reaction.Update_Conversion(*ETA, *T, Delta_t);
    }
}

std::pair<std::vector<long double>, std::vector<long double>> Combustion_Problem::Iterate()
{
    Setup_Solver();
    Solve_and_Update_State();

    std::vector<long double> Temperature(T_VECTOR);
    std::vector<long double> Conversion(ETA_VECTOR);

    return std::pair<std::vector<long double>, std::vector<long double>> (Temperature, Conversion);
}

void Combustion_Problem::Set_Initial_Conditions(
    std::vector<long double> T_Array,
    std::vector<long double> Eta_Array
)
{
    std::copy(T_Array.begin(), T_Array.end(), T_VECTOR.begin());
    std::copy(Eta_Array.begin(), Eta_Array.end(), ETA_VECTOR.begin());
}

bool Combustion_Problem::Combustion_Not_Completed()
{
    bool flag = true;

    for (ETA = ETA_VECTOR.begin(); ETA < (ETA_VECTOR.end()-1) && flag; ETA++)

        flag = (*ETA) > 1.0l - 0.00001l;

    return !flag;
}

void Combustion_Problem::Write_to_File(std::ofstream &file)
{
    file << "Number of Grid Points :\t" << N << std::endl;

    file << std::endl;

    file << "Delta x :\t" << Delta_x << "\t m" << std::endl;
    file << "Delta t :\t" << Delta_t << "\t s" << std::endl; 

    file << std::endl;

    file << "Ignition Temperature :\t" << Ignition_Temperature << "\t K" << std::endl;

    file << std::endl;
}

