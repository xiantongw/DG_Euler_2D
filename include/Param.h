#ifndef PARAM_H
#define PARAM_H
// Simulation Parameter is stored in a Struct
typedef struct Param{
    double cfl;
    double mach_inf;
    double attack_angle;
    double gamma;
    double p_inf;
    double T_inf;
    double R;
    double eps;
    int dnOutput;
    int MAXITER;
    std::string bound0;
    std::string bound1;
    std::string bound2;
    std::string bound3;
    std::string mesh_file;
    int order;
    int order_geo;
} Param;

#endif
