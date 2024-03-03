#ifndef _INITIAL_CONDITIONS_H_
#define _INITIAL_CONDITIONS_H_

#include <math.h>
#include "vectors.h"
#include "particles.h"
#include "log.h"
#include "psrand.h"

static inline void initialize_random(particles parts, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax){
    LOG("initializing particles randomly in box %g <= x <= %g, %g <= y <= %g, %g <= z <= %g",xmin,xmax,ymin,ymax,zmin,zmax);
    vec3* pos = parts.pos;
    vec3* vel = parts.vel;
    int n_particles = parts.n_particles;

    for (int i = 0; i < n_particles; i++){
        vec3 p_i = v3(frand_range(xmin,xmax),
                        frand_range(ymin,ymax),
                        frand_range(zmin,zmax));
        vec3 v_i = v3zero();
        pos[i] = p_i;
        vel[i] = v_i;
    }
    LOG("initialized particles randomly");
}

static inline void initialize_random_sphere(particles parts, float r){
    LOG("initializing particles randomly in sphere radius = %g",r);
    vec3* pos = parts.pos;
    vec3* vel = parts.vel;
    int n_particles = parts.n_particles;

    for (int i = 0; i < n_particles; i++){
        vec3 p_i = v3(frand_range(-r,r),
                        frand_range(-r,r),
                        frand_range(-r,r));
        while (v3len(p_i) > r){
            p_i = v3(frand_range(-r,r),
                        frand_range(-r,r),
                        frand_range(-r,r));
        }
        vec3 v_i = v3zero();
        pos[i] = p_i;
        vel[i] = v_i;
    }
    LOG("initialized particles randomly in sphere");
}

#endif