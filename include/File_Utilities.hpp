#ifndef __FILE_UTIL__
#define __FILE_UTIL__

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <string.h>

#include "Combustion_Problem.hpp"

extern Substance Preheat_Zone_Aluminium;
extern Substance Preheat_Zone_Nickel;
extern Substance Preheat_Zone_Argon;

extern Substance Post_Combustion_Zone_NiAl;
extern Substance Post_Combustion_Zone_Argon;

extern Substance Reaction_Zone_Aluminium;
extern Substance Reaction_Zone_Nickel;
extern Substance Reaction_Zone_Argon;
extern Substance Reaction_Zone_NiAl;

extern Coated_Particle Pre_Heat_Zone_Ni_Coated_Al_Particle;
extern Coated_Particle Reaction_Zone_Ni_Coated_Al_Particle;
extern Coated_Particle Reaction_Zone_NiAl_Particle;
extern Coated_Particle Post_Combustion_Zone_NiAl_Particle;

void Make_Solution_Folder();

void Save_Conversion_Data(std::vector<std::vector<long double>> ETA_MATRIX);
void Save_Temperature_Data(std::vector<std::vector<long double>> T_MATRIX);

void Save_Combustion_Config(
    unsigned int Number_of_Time_Steps,
    long double Particle_Volume_Fraction,
    Combustion_Pellet Pellet,
    Combustion_Problem Problem
);

#endif