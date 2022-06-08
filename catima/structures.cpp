#include "structures.h"
#include "catima/nucdata.h"
#include <algorithm>
#include <cmath>


namespace catima{

bool operator==(const Projectile &a, const Projectile&b){
    if( (a.A==b.A) && (a.Z==b.Z) && (a.Q==b.Q)){
        return true;
        }
    else
        return false;
}

bool operator==(const Material &a, const Material&b){
    if(std::fabs(a.density() - b.density())> 1e-6)return false;
    if(a.ncomponents() != b.ncomponents())return false;
    if(a.I() != b.I())return false;
    for(int i=0;i<a.ncomponents();i++){
        if(a.atoms[i].stn != b.atoms[i].stn)return false;
        if(a.atoms[i].A != b.atoms[i].A)return false;
        if(a.atoms[i].Z != b.atoms[i].Z)return false;
    }
    if(a.molar_mass != b.molar_mass)return false;
    return true;
}


Material::Material(std::initializer_list<std::array<double,3>>list,double _density, double _ipot, double mass):rho(_density),i_potential(_ipot){
    std::initializer_list<std::array<double,3>>::iterator it;
    atoms.reserve(list.size());
    for ( it=list.begin(); it!=list.end(); ++it){
        add_element((*it)[0],(*it)[1],(*it)[2]);
    }

    if(mass!=0.0){
        molar_mass=mass;
        }
    else{
        calculate(); // calculate if needed, ie average molar mass
        }
}

Material::Material(double _a, int _z, double _rho, double _th, double _ipot):rho(_rho),th(_th),i_potential(_ipot){
    add_element(_a,_z,1.0);
}

void Material::add_element(double _a, int _z, double _stn){
    double a = (_a>0)?_a:element_atomic_weight(_z);
    atoms.push_back({a,_z,_stn});
    molar_mass += _stn*a;
}

void Material::calculate(){
    if(std::all_of(atoms.begin(),atoms.end(),[](const Target &t){return t.stn<1.0;})){
        double sum = 0;
        for(auto& e: atoms){
            sum+= e.stn/e.A;
        }
        molar_mass = 1.0/sum;
    }
}

void Layers::add(Material m){
    materials.push_back(m);
}

void Layers::add(const Layers& l){
    for(auto m: l.get_materials()){
      add(m);
    }
}

double Layers::thickness() const {
    double sum = 0;
    for(auto &m : materials){
        sum += m.thickness();
    }
    return sum;
}

double Layers::thickness_cm() const {
    double sum = 0;
    for(auto &m : materials){
        sum += m.thickness_cm();
    }
    return sum;
}

Layers operator+(const Layers &a, const Layers&b){
    Layers res;
    for(auto &e:a.materials){
        res.add(e);
    }

    for(auto &e:b.materials){
        res.add(e);
    }
    return res;
}

Layers operator+(const Layers &a, const Material &m){
    Layers res;
    for(auto &e:a.materials){
        res.add(e);
    }
    res.add(m);
    return res;
}


}
