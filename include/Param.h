#ifndef PARAM_H
#define PARAM_H
// Simulation Parameter is stored in a Struct
typedef struct Param{
    double cfl;
    double mach_inf;
    double attack_angle;
    double gamma;
    double p_inf;
    char *bound0;
    char *bound1;
    char *bound2;
    char *bound3;
} Param;

#endif
