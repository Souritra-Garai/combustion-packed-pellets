#ifndef __SUBSTANCE__
#define __SUBSTANCE__

#include <fstream>

class Substance
{
    private:
        
        long double Density;
        long double Heat_Capacity;
        long double Heat_Conductivity;

        long double Molecular_Weight;

    public:
        
        Substance(
            long double Density,
            long double Heat_Capacity,
            long double Heat_Conductivity,
            long double Molecular_Weight
        );
        
        long double Get_Density();
        long double Get_Heat_Capacity();
        long double Get_Heat_Conductivity();
        
        long double Get_Molecular_Weight();

        void Write_to_File(std::ofstream &File_Handler, const char *Name);
};

#endif