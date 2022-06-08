/*
 *  Author: Andrej Prochazka
 *  Copyright(C) 2017
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/// \file calculations.h
#ifndef CALCULATIONS_H
#define CALCULATIONS_H
#include <complex>
#include "catima/structures.h"
#include "catima/config.h"

namespace catima{

    /**
      * returns nuclear stopping power for projectile-target combination
      */
    double dedx_n(const Projectile &p, const Target &t);
    double dedx_n(const Projectile &p, const Material &mat); 
    
    /**
      * returns energy loss straggling
      */
    double dedx_variance(const Projectile &p, const Target &t, const Config &c=default_config);

    /**
      * returns reduced energy loss unit for projectile-target combination
      */
    double reduced_energy_loss_unit(const Projectile &p, const Target &t);
    

    /**
     * @brief bethek_dedx_e - electronics stopping power
     * @return stopping power
     */
    double bethek_dedx_e(const Projectile &p,const Target &t, const Config &c=default_config, double I=0.0);
    double bethek_dedx_e(const Projectile &p,const Material &mat, const Config &c=default_config);

    /** 
      * calculates barkas effect
      */
    double bethek_barkas(double zp_eff,double eta, double zt);
    
    /** 
      * calculates density effect
      */
    double bethek_density_effect(double beta, int zt);
    
    /**
      * calculates lindhard correction for energy loss calculation
      */
    double bethek_lindhard(const Projectile &p);

    /**
      * calculates lindhard correction for energy loss straggling calculation
      */
    double bethek_lindhard_X(const Projectile &p);

    /**
      * calculates pair production stopping power
      */
    double pair_production(const Projectile &p, const Target &t);

    /**
      * calculates bremsstrahlung stopping power
      */
    double bremsstrahlung(const Projectile &p, const Target &t);

    /**
      * returns linhard correction (L) calulated from tabulated precalculated data
      * if energy is less than minimal calculated energy the LS coefficient of at  minimal 
      * calculated energy is returned and similar for highest caclulated energy limit
      */
    double precalculated_lindhard(const Projectile &p);

    /**
      * returns linhard energy loss straggling correction (X) calulated from tabulated precalculated data
      * if energy is less than minimal calculated energy the X coefficient of at  minimal 
      * calculated energy is returned and similar for highest caclulated energy limit
      */
    double precalculated_lindhard_X(const Projectile &p);

    /**
      * this function is not used and is not tested 
      */
    double energy_straggling_firsov(double z1,double energy, double z2, double m2);
     
    /**
      * electronic energy loss for low energy, should be like SRIM
      */ 
    double sezi_dedx_e(const Projectile &p, const Material &mat, const Config &c=default_config);

    //constexpr double Es2_FR =2*PI/fine_structure* electron_mass * electron_mass;
    constexpr double Es2_FR = 198.81;
    
    /**
     * angular scattering power in form of da^2/dx in units rad^2/ g/cm^2 
     * @param p - Projectile class
     * @param t - Target class
     * @Es2 - energy constant squared, default is 14.1^2 = 198.81
     */
    double angular_scattering_power(const Projectile &p, const Target &t, double Es2=Es2_FR);
    
    /**
     * angular scattering power in form of da^2/dx in units rad^2/ g/cm^2 
     * @param p - Projectile class
     * @param t - Material class
     * @Es2 - energy constant squared, default is 14.1^2 = 198.81
     */
    double angular_scattering_power(const Projectile &p, const Material &material, double Es2=Es2_FR);
    double angular_scattering_power_xs(const Projectile &p, const Material &mat, double p1, double beta1, double Es2=225.0);
    /**
      * returns radiation length of the (M,Z) material
      * for certain z the radiation length is tabulated, otherwise calculated
      * @param z - proton number of target
      * @param m - weight of the target
      * @return radiation length in g/cm^2
      */
    double radiation_length(int z, double m);

    /**
      * returns radiation length of the Material class
      * radiation length if calculated if not specified in Material
      * or return specified radiation length      
      * @param Material - Material class
      * @return radiation length in g/cm^2
      */
    double radiation_length(const Material &material);

    /**
      * returns radiation length of the Material class
      * radiation length if calculated if not specified in Material
      * or return specified radiation length      
      * @param Material - Material class
      * @return radiation length in g/cm^2
      */
    double scattering_length(const Material &material);
    
    /** returns effective Z of the projectile
      * @param p - Projectile class
      * @param t - Target class
      * @param c - Configuration, the z effective will be calculated according to c.z_effective value
      * @return - z effective
      */ 
    double z_effective(const Projectile &p, const Target &t, const Config &c=default_config);

    /**
      * calculates effective charge
      * @param z - proton number of projectile
      * @param beta - velocity of projectile
      * @return effective charge
      */
    double z_eff_Pierce_Blann(double z, double beta);

    /**
      * calculates effective charge
      * @param pz - proton number of projectile
      * @param beta - velocity of projectile
      * @param tz - proton number of target material
      * @return effective charge
      */
    double z_eff_Anthony_Landford(double pz, double beta, double tz);

    /**
      * calculates effective charge
      * @param pz - proton number of projectile
      * @param beta - velocity of projectile
      * @param tz - proton number of target material
      * @return effective charge
      */
    double z_eff_Hubert(double pz, double E, double tz);
    
    /**
      * calculates effective charge
      * @param pz - proton number of projectile
      * @param beta - velocity of projectile
      * @param tz - proton number of target material
      * @return effective charge
      */
    double z_eff_Winger(double pz, double beta, double tz);

    /**
      * calculates effective charge
      * @param pz - proton number of projectile
      * @param beta - velocity of projectile
      * @param tz - proton number of target material
      * @return effective charge
      */
    double z_eff_global(double pz, double E, double tz);

    /**
      * calculates effective charge
      * @param pz - proton number of projectile
      * @param beta - velocity of projectile
      * @param tz - proton number of target material
      * @return effective charge
      */
    double z_eff_Schiwietz(double pz, double beta, double tz);

    /**
      * calculates effective charge
      * @param pz - proton number of projectile
      * @param beta - velocity of projectile
      * @param tz - proton number of target material
      * @return effective charge
      */
    double z_eff_atima14(double pz, double beta, double tz);

    

    //helper
    double gamma_from_T(double T);
    double beta_from_T(double T);
    double p_from_T(double T, double M);
    std::complex<double> lngamma( const std::complex<double> &z );
    std::complex<double> hyperg(const std::complex<double> &a,
                                                const std::complex<double> &b,
                                                const std::complex<double> &z);

    inline double power(double x, double y){
        return exp(log(x)*y);
    }

    inline double ipow(double x, int y){
      if(y==0)return 1.0;
      else if(y==1) return x;      
      else if(y==2) return x*x;      
      else return std::pow(x,y);
    }

}
#endif
