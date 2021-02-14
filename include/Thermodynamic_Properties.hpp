#ifndef __THERMO_PROPS__
#define __THERMO_PROPS__

#include <math.h>

#define lambda_Al_P 220.0
#define lambda_Al_R 130.0

#define lambda_Ni_P 66.0
#define lambda_Ni_R 80.0

#define lambda_Ar_P 0.016
#define lambda_Ar_R 0.055

#define lambda_p_P 137.0
#define lambda_p_R 108.0

#define lambda_NiAl 115.0

#define rho_Al      2700.0
#define rho_Ni      8908.0
#define rho_NiAl    5900.0
#define rho_p       5430.0
#define rho_Ar_P    0.89
#define rho_Ar_R    0.37

#define Cp_Al_P 1060.0
#define Cp_Al_R 1176.0

#define Cp_Ni_P 440.0
#define Cp_Ni_R 670.0

#define Cp_Ar_P 520.0
#define Cp_Ar_R 520.0

#define Cp_p_P 613.0
#define Cp_p_R 760.0

#define Cp_NiAl 717.0

float get_lambda_m(float lambda_p, float lambda_f, float phi_p, float alpha_p_1);

float get_Cp_m(float phi_p, float Cp_p, float Cp_f, float rho_f);

float get_rho_m(float phi_p, float rho_f);

#endif