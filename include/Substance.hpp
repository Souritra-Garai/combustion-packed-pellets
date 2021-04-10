#ifndef __SUBSTANCE__
#define __SUBSTANCE__

class Substance
{
    private:
        
        long double Density;
        long double Heat_Capacity;
        long double Heat_Conductivity;

    public:
        
        Substance(long double DENSITY, long double HEAT_CAPACITY, long double HEAT_CONDUCTIVITY);
        
        long double Get_Density();
        long double Get_Heat_Capacity();
        long double Get_Heat_Conductivity();
};

#endif