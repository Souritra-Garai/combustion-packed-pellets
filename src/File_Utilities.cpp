#include "File_Utilities.hpp"

time_t rawtime;
struct tm * timeinfo;
char folder[100] = "solutions/";
char T_file[100], Eta_file[100], config_file[100];

void Make_Solution_Folder()
{
    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    char *curr_time = asctime(timeinfo);
    curr_time[strlen(curr_time) - 1] = '\0';

    strcat(folder, asctime(timeinfo));

    strcpy(T_file, folder);
    strcpy(Eta_file, folder);
    strcpy(config_file, folder);

    mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    strcat(T_file, "/Temperature.csv");
    strcat(Eta_file, "/Conversion.csv");
    strcat(config_file, "/Combustion_Config.txt");
}

void Save_Temperature_Data(std::vector<std::vector<long double>> T_MATRIX)
{
    std::ofstream file(T_file);

    for (auto T_VECTOR = T_MATRIX.begin(); T_VECTOR < T_MATRIX.end(); T_VECTOR++)
    {
        for (auto T = T_VECTOR->begin(); T < T_VECTOR->end(); T++)
        {
            file << std::to_string((*T)) << ',';
        }

        file << std::endl;
    }

    file.close();
}

void Save_Conversion_Data(std::vector<std::vector<long double>> T_MATRIX)
{
    std::ofstream file(Eta_file);

    for (auto T_VECTOR = T_MATRIX.begin(); T_VECTOR < T_MATRIX.end(); T_VECTOR++)
    {
        for (auto T = T_VECTOR->begin(); T < T_VECTOR->end(); T++)
        {
            file << std::to_string((*T)) << ',';
        }

        file << std::endl;
    }

    file.close();
}

void Save_Combustion_Config(
    unsigned int n,
    long double Particle_Volume_Fraction,
    Combustion_Pellet Pellet,
    Combustion_Problem my_Problem
)
{
    std::ofstream txt_file(config_file);

    txt_file << "Number of Time Steps :\t" << n+1 << std::endl;

    my_Problem.Write_to_File(txt_file);

    Pellet.Write_to_File(txt_file, "Ni Coated Al Pellet degassed with Ar");

    txt_file << "Pellet Volume Fraction Occupied by Solid Particles :\t" << Particle_Volume_Fraction << std::endl;

    txt_file << std::endl;

    Pre_Heat_Zone_Ni_Coated_Al_Particle.Write_to_File(txt_file, "Pre-Heat Zone Ni Coated Al Particle");
    Reaction_Zone_Ni_Coated_Al_Particle.Write_to_File(txt_file, "Reaction Zone Ni Coated Al Particle");
    Reaction_Zone_NiAl_Particle.Write_to_File(txt_file, "Reaction Zone Product NiAl Particle");
    Post_Combustion_Zone_NiAl_Particle.Write_to_File(txt_file, "Post-Combustion Zone Product NiAl Particle");

    Preheat_Zone_Aluminium.Write_to_File(txt_file, "Pre-Heat Zone Al");
    Preheat_Zone_Nickel.Write_to_File(txt_file, "Pre-Heat Zone Ni");
    Preheat_Zone_Argon.Write_to_File(txt_file, "Pre-Heat Zone Ar");

    Reaction_Zone_Aluminium.Write_to_File(txt_file, "Reaction Zone Al");
    Reaction_Zone_Nickel.Write_to_File(txt_file, "Reaction Zone Ni");
    Reaction_Zone_NiAl.Write_to_File(txt_file, "Reaction Zone NiAl");
    Reaction_Zone_Argon.Write_to_File(txt_file, "Reaction Zone Ar");

    Post_Combustion_Zone_NiAl.Write_to_File(txt_file, "Post-Combustion Zone NiAl");
    Post_Combustion_Zone_Argon.Write_to_File(txt_file, "Post-Combustion Zone Ar");

    txt_file.close();
}