#ifndef __THERMO_PROPS__
#define __THERMO_PROPS__

#include <math.h>

// Preheat Zone
#define lambda_p_P  137.0
#define lambda_Ar_P 0.016

#define rho_p_P     5430.0
#define rho_Ar_P    0.89

#define Cp_p_P  613.0
#define Cp_Ar_P 520.0

// Post Combustion Zone
#define lambda_Ar_PC    0.055
#define lambda_NiAl_PC  115.0

#define rho_Ar_PC   0.37
#define rho_NiAl_PC 5900.0

#define Cp_Ar_PC    520.0
#define Cp_NiAl_PC  717.0

// Reaction Zone
#define lambda_p_R      108.0
#define lambda_NiAl_R   115.0
#define lambda_Ar_R     0.055

#define rho_p_R     5430.0
#define rho_NiAl_R  5900.0
#define rho_Ar_R    0.89

#define Cp_p_R      760.0
#define Cp_NiAl_R   717.0
#define Cp_Ar_R     520.0

// Calculates thermal conductivity of packed particles
// and fluid mixture (Ref. Section 3.1.1)
// phi_p - volume fraction of packed particles
// alpha_p_ME - volume fraction of packed particles in Maxwell - Eucken structure
// lambda_p - thermal conductivity of particles
// lambda_f - thermal conductivity of fluid
float calc_lambda_m(float phi_p, float alpha_p_ME, float lambda_p, float lambda_f);

// Calculates specific heat capacity of packed particles
// and fluid mixture (Ref. Section 3.3.1)
// phi_p - volume fraction of packed particles
// rho_p - density of particles
// rho_f - density of fluid
// Cp_p - specific heat capacity of particles
// Cp_f - specific heat capacity of fluid
float calc_Cp_m(float phi_p, float rho_p, float rho_f, float Cp_p, float Cp_f);

// Calculates density of packed particles and fluid
// mixture (Ref. Section 3.2.1)
// phi_p - volume fraction of packed particles
// rho_p - density of packed particles
// rho_f - density of fluid
float calc_rho_m(float phi_p, float rho_p, float rho_f);

#endif