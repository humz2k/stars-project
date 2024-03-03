#ifndef _SMOOTHING_KERNEL_H_
#define _SMOOTHING_KERNEL_H_

#include <math.h>
#include "vectors.h"

#define SQRT_PI 1.77245385091
#define PI_3_2 5.56832799683

/*For SPH, the `smoothed` refers to the fact that we
distribute the mass of a particle smoothly around the
particle using some function.

For this code, we will use a Gaussian smoothing kernel:

    W(vec(r), h) = [1 / (h^3 * pi^(3/2))] * exp(-|r|^2 / h^2)

where |r| is the distance from the particle, and h is the
smoothing distance.*/
static inline float W(vec3 r, float h){
    float r2 = v3len2(r);
    float a = (1.0f / (h * SQRT_PI));
    float b = expf(-r2 / (h * h));
    float w = a * a * a * b;
    return w;
}

/*We also need the gradient of W:
Grad(W(vec(r), h)) = [-2/(h^5 * pi^(3/2))] * exp(-|r|^2 / h^2) * vec(r)*/
static inline vec3 gradW(vec3 r, float h){
    float r2 = v3len2(r);
    float h2 = h * h;
    float h5 = h2 * h2 * h;
    float a = -(2.0f / (h5 * PI_3_2));
    float b = expf(-r2 / (h * h));
    vec3 gw = v3fmul(r,a * b);
    return gw;
}

#endif