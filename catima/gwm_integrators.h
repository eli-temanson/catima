#ifndef GWM_INTEGRATORS_H
#define GWM_INTEGRATORS_H

#include "catima/catima.h"

namespace catima {

    double integrate_energyloss(Projectile& proj, const Material& mat, const Config& c=default_config);

    double reverse_integrate_energyloss(Projectile& proj, const Material& mat, const Config& c=default_config);
}

#endif