#include "Thermo_Physical_Properties/Substance.hpp"

Substance::Substance(long double rho, long double C, long double lambda)
{
    Density = rho;
    Heat_Capacity = C;
    Heat_Conductivity = lambda;
}

long double Substance::Get_Density()            {return Density;}
long double Substance::Get_Heat_Capacity()      {return Heat_Capacity;}
long double Substance::Get_Heat_Conductivity()  {return Heat_Conductivity;}
