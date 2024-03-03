#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "vectors.h"
#include "log.h"

/*holds particle positions, velocities, accelerations, densities and mass*/
typedef struct particles{
    int n_particles;
    float mass;
    vec3* pos;
    vec3* vel;
    vec3* acc;
    float* P_rho2;
} particles;

/*allocates n particles*/
static inline particles particles_alloc(int n_particles, float mass){
    LOG("allocating %d particles with mass %g",n_particles,mass);

    size_t vec3size = sizeof(vec3) * n_particles;
    size_t floatsize = sizeof(float) * n_particles;
    size_t total = 0;

    particles out;
    out.mass = mass;
    out.n_particles = n_particles;
    assert(out.pos = (vec3*)malloc(vec3size));
    memset(out.pos,0,vec3size);
    LOG("allocated %lu bytes for pos",vec3size);
    total += vec3size;

    assert(out.vel = (vec3*)malloc(vec3size));
    memset(out.vel,0,vec3size);
    LOG("allocated %lu bytes for vel",vec3size);
    total += vec3size;

    assert(out.acc = (vec3*)malloc(vec3size));
    memset(out.acc,0,vec3size);
    LOG("allocated %lu bytes for acc",vec3size);
    total += vec3size;

    assert(out.P_rho2 = (float*)malloc(floatsize));
    memset(out.pos,0,floatsize);
    LOG("allocated %lu bytes for P_rho2",floatsize);
    total += floatsize;

    LOG("allocated %lu bytes total",total);

    LOG("allocated particles");
    return out;
}

/*deletes allocated particles*/
static inline void particles_destroy(particles p){
    LOG("destroying particles");
    free(p.pos);
    free(p.vel);
    free(p.acc);
    free(p.P_rho2);
    LOG("destroyed particles");
}

static inline void particles_debug(particles p){
    for (int i = 0; i < p.n_particles; i++){
        printf("%.2f %.2f %.2f | %.2f %.2f %.2f | %.2f %.2f %.2f | %.2f |\n",p.pos[i].x,p.pos[i].y,p.pos[i].z,p.vel[i].x,p.vel[i].y,p.vel[i].z,p.acc[i].x,p.acc[i].y,p.acc[i].z,p.P_rho2[i]);
    }
}

#endif