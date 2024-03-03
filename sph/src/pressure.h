#ifndef _DENSITY_H_
#define _DENSITY_H_

#include <math.h>
#include "vectors.h"
#include "smoothing_kernel.h"
#include "particles.h"

/*Pressure according to polytrope equation of state:
    `P = k * rho^(1+1/n)`
*/
static inline float getPressure(float rho,  // density
                                float k,    // equation of state constant
                                float n     // polytropic index
                                ){
    float P = k * powf(rho,1.0f + (1.0f/n));
    return P;
}

#endif