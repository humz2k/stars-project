#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "vectors.h"
#include "particles.h"
#include "acceleration.h"

/*updates velocites based on cached acceleration*/
static void updateVelocities(particles parts, float dt){
    int n_particles = parts.n_particles;
    vec3* vel = parts.vel;
    vec3* acc = parts.acc;

    #pragma omp parallel for
    for (int i = 0; i < n_particles; i++){
        vec3 v_i = vel[i];
        vec3 a_i = acc[i];
        v_i = v3v3add(v_i,v3fmul(a_i,dt));
        vel[i] = v_i;
    }
}

/*updates positions based on cached velocities*/
static void updatePositions(particles parts, float dt){
    int n_particles = parts.n_particles;
    vec3* vel = parts.vel;
    vec3* pos = parts.pos;

    #pragma omp parallel for
    for (int i = 0; i < n_particles; i++){
        vec3 v_i = vel[i];
        vec3 p_i = pos[i];
        p_i = v3v3add(p_i,v3fmul(v_i,dt));
        pos[i] = p_i;
    }
}

/*does a leap frog integration step*/
static void leapfrog_timestep(particles parts, float dt, float h, float k, float n, float lambda, float nu){
    updateVelocities(parts,dt * 0.5f);
    updatePositions(parts,dt);
    cacheP_rho2(parts,h,k,n);
    cacheAcc(parts,h,k,n,lambda,nu);
    updateVelocities(parts,dt * 0.5f);
}

/*initializes the leapfrog timestepper*/
static void leapfrog_init(particles parts, float h, float k, float n, float lambda, float nu){
    LOG("initializing leapfrog");
    cacheP_rho2(parts,h,k,n);
    cacheAcc(parts,h,k,n,lambda,nu);
}

#endif