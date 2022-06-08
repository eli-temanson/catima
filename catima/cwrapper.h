#ifndef CATIMA_CWRAPPER
#define CATIMA_CWRAPPER

#ifdef __cplusplus
extern "C" {
#endif

struct CatimaResult{
        double Ein;
        double Eout;
        double Eloss;
        double range;
        double dEdxi;
        double dEdxo;
        double sigma_E;
        double sigma_a;
        double sigma_r;
        double tof;
};


    enum z_eff_type {
        none = 0,
        pierce_blann = 1,
        anthony_landorf = 2,
        hubert = 3,
        winger = 4,
        schiwietz = 5,
        global = 6,
        atima14 = 7
    };

struct CatimaConfig {
        char z_effective;
};

extern struct CatimaConfig catima_defaults;

typedef struct CatimaResult CatimaResult;

CatimaResult catima_calculate(double pa, int pz, double T, double ta, double tz, double thickness, double density);
double catima_Eout(double pa, int pz, double T, double ta, double tz, double thickness, double density);
double catima_range(double pa, int pz, double T, double ta, double tz);
double catima_range_straggling(double pa, int pz, double T, double ta, double tz);
double catima_angular_straggling_from_E(double pa, int pz, double Tin, double Tout,double ta, double tz);
double catima_energy_straggling_from_E(double pa, int pz, double Tin, double Tout,double ta, double tz);
double atomic_weight(int i);
double catima_nonreaction_rate(double pa, int pz, double T, double ta, double tz, double thickness);

#ifdef __cplusplus
}
#endif

#endif
