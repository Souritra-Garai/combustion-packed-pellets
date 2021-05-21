#include "Combustion_Problem.hpp"
// #include <iostream>

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
) : Pellet(P), Combustion_Reaction(R)
{
    N = n;

    T_VECTOR = new long double[N];
    
    ETA_VECTOR = new long double[N];
    
    Delta_t = Dt;

    Ignition_Temperature = T_ign;
    Adiabatic_Combustion_Temperature = T_ad;

    Delta_x = Pellet.Get_Pellet_Length() / (N - 1);
}

bool Combustion_Problem::In_Post_Combustion_Zone(long double Conversion, long double Temperature)
{
    return (Conversion > 1.0l - 0.00001l);
    // return (Temperature >= Adiabatic_Combustion_Temperature) && In_Range_Equation_Iterators();
}

bool Combustion_Problem::In_Reaction_Zone(long double Temperature)
{
    return (Temperature >= Ignition_Temperature);
}

void Combustion_Problem::Solve_and_Update_State()
{   
    long double T_VECTOR_NXT[N];

    // T_VECTOR_NXT[0] = Adiabatic_Combustion_Temperature; // Isothermal BC
    // T_VECTOR_NXT[0] = T_VECTOR[1] - Delta_x * 0.0l / (Pellet.Get_Cross_Section_Area() * Pellet.Get_Post_Combustion_Zone_Diffusion_Term_Coeff());  // Adiabatic BC 
    T_VECTOR_NXT[0] = T_VECTOR[1] - Delta_x * (Pellet.Get_Flat_Surface_Heat_Loss_Term(T_VECTOR[0])) / (Pellet.Get_Cross_Section_Area() * Pellet.Get_Post_Combustion_Zone_Diffusion_Term_Coeff()); // Heat Loss BC

    #pragma omp parallel for

        for (int i = 1; i < N-1; i++)
        {
            if (T_VECTOR[i] >= Ignition_Temperature)
            {
                if (ETA_VECTOR[i] >= 1.0l - 0.00001l)
                {
                    T_VECTOR_NXT[i] = T_VECTOR[i] + (Delta_t / Pellet.Get_Post_Combustion_Zone_Transient_Term_Coeff()) * (
                        Pellet.Get_Post_Combustion_Zone_Diffusion_Term_Coeff() * (T_VECTOR[i+1] - 2 * T_VECTOR[i] + T_VECTOR[i-1]) / (Delta_x * Delta_x) -
                        Pellet.Get_Lateral_Surface_Heat_Loss_Term(T_VECTOR[i])
                    );
                }

                else
                {
                    long double omega = Combustion_Reaction.Get_Omega(T_VECTOR[i], ETA_VECTOR[i]);

                    T_VECTOR_NXT[i] = T_VECTOR[i] + (Delta_t / Pellet.Get_Reaction_Zone_Transient_Term_Coeff()) * (
                        Pellet.Get_Reaction_Zone_Diffusion_Term_Coeff() * (T_VECTOR[i+1] - 2 * T_VECTOR[i] + T_VECTOR[i-1]) / (Delta_x * Delta_x) +
                        Combustion_Reaction.Get_Energy_Gen_Rate(omega) -
                        Pellet.Get_Lateral_Surface_Heat_Loss_Term(T_VECTOR[i])
                    );

                    ETA_VECTOR[i] += omega * Delta_t;
                }
            }

            else
            {
                T_VECTOR_NXT[i] = T_VECTOR[i] + (Delta_t / Pellet.Get_Pre_Heat_Zone_Transient_Term_Coeff()) * (
                    Pellet.Get_Pre_Heat_Zone_Diffusion_Term_Coeff() * (T_VECTOR[i+1] - 2 * T_VECTOR[i] + T_VECTOR[i-1]) / (Delta_x * Delta_x) -
                    Pellet.Get_Lateral_Surface_Heat_Loss_Term(T_VECTOR[i])
                );
            }   
        }

    // T_VECTOR_NXT[N-1] = T_VECTOR[N-2] - Delta_x * 0.0l / (Pellet.Get_Cross_Section_Area() * Pellet.Get_Post_Combustion_Zone_Diffusion_Term_Coeff());  // Adiabatic BC 
    T_VECTOR_NXT[N-1] = T_VECTOR[N-2] - Delta_x * (Pellet.Get_Flat_Surface_Heat_Loss_Term(T_VECTOR[N-1])) / (Pellet.Get_Cross_Section_Area() * Pellet.Get_Post_Combustion_Zone_Diffusion_Term_Coeff()); // Heat Loss BC
    
    #pragma omp parallel for
    
        for (int i = 0; i < N; i++) T_VECTOR[i] = T_VECTOR_NXT[i];
}

std::pair<std::vector<long double>, std::vector<long double>> Combustion_Problem::Iterate()
{
    Solve_and_Update_State();

    std::vector<long double> Temperature(N);
    std::vector<long double> Conversion(N);

    for (int i = 0; i < N; i++)
    {
        Temperature[i] = T_VECTOR[i];
        Conversion[i]  = ETA_VECTOR[i];
    }

    return std::pair<std::vector<long double>, std::vector<long double>> (Temperature, Conversion);
}

void Combustion_Problem::Set_Initial_Conditions(
    std::vector<long double> T_Array,
    std::vector<long double> Eta_Array
)
{
    for (int i = 0; i < N; i++)
    {
        T_VECTOR[i] = T_Array[i];
        ETA_VECTOR[i] = Eta_Array[i];
    }
}

bool Combustion_Problem::Combustion_Not_Completed()
{
    bool flag = true;

    for (int i = 0; i < N && flag; i++)

        flag = ETA_VECTOR[i] > 1.0l - 0.00001l;

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

Combustion_Problem::~Combustion_Problem()
{
    delete[] T_VECTOR;
    delete[] ETA_VECTOR;

    // printf("\nCombustion Probelem destructor was called\n");
}