#ifndef _ACCELERATION_H_
#define _ACCELERATION_H_

#include <math.h>
#include "vectors.h"
#include "particles.h"
#include "pressure.h"
#include "smoothing_kernel.h"

/*Gets the density at position `p_i`*/
static inline float getRho(vec3 p_i, particles parts, float h){
    float mass = parts.mass;
    vec3* pos = parts.pos;
    int n_particles = parts.n_particles;
    float rho = 0.0f;
    for (int j = 0; j < n_particles; j++){
        vec3 p_j = pos[j];
        float w = W(v3v3sub(p_i,p_j),h);
        rho += w * mass;
    }
    return rho;
}

/*caches `P_i/rho_i^2` of all particles*/
static inline void cacheP_rho2(particles parts, float h, float k, float n){
    // read pointers/data from `parts`
    int n_particles = parts.n_particles;
    vec3* pos = parts.pos;
    float* P_rho2 = parts.P_rho2;

    // loop through particles `i`
    #pragma omp parallel for
    for (int i = 0; i < n_particles; i++){
        // load `p_i`
        vec3 p_i = pos[i];

        // calculate `rho_i` and `P_i`
        float rho = getRho(p_i,parts,h);
        float P = getPressure(rho,k,n);

        // cache `P_i/rho_i^2`
        P_rho2[i] = P / (rho * rho);
    }
}
/*Calculates force from external potential:
    `vec(f) = -lambda * vec(r) - nu * vec(v)`
*/
static inline vec3 calcExternalForce(vec3 p_i, vec3 v_i, float lambda, float nu){
    return v3v3add(v3fmul(p_i,-lambda),v3fmul(v_i,-nu));
}

/*calculates acceleration at particle p_i:
    `dv_i/dt = -sum_{j,j != i}[ m * ((P_i / rho_i^2) + (P_j / rho_j^2)) * grad(W(r_i - r_j, h)) + vec(f_i) ]`
*/
static inline vec3 calcAcc(int i, particles parts, float h, float k, float n, float lambda, float nu){

    assert(i < parts.n_particles);

    // read pointers/data from `parts`
    vec3* pos = parts.pos;
    vec3* acc = parts.acc;
    vec3* vel = parts.vel;
    float* P_rho2 = parts.P_rho2;
    float mass = parts.mass;
    int n_particles = parts.n_particles;

    // read pos/rho for `i`
    vec3 p_i = pos[i];
    vec3 v_i = vel[i];
    float P_rho2_i = P_rho2[i];

    // initialize acceleration at `0`
    vec3 a_i = v3zero();

    // loop through all particles `j, j != i`
    for (int j = 0; j < n_particles; j++){
        if (i == j)continue;

        // grab `p_j` and `rho_j`
        vec3 p_j = pos[j];
        float P_rho2_j = P_rho2[j];

        // calculate `gradw`
        vec3 gradw = gradW(v3v3sub(p_i,p_j),h);

        // calculate acceleration magnitude
        float mag_a_i = -mass * (P_rho2_i + P_rho2_j);

        a_i = v3v3add(a_i,v3fmul(gradw,mag_a_i));
    }

    return v3v3add(a_i,calcExternalForce(p_i,v_i,lambda,nu));
}

/*caches acceleration for all particles*/
static inline void cacheAcc(particles parts, float h, float k, float n, float lambda, float nu){
    vec3* acc = parts.acc;
    int n_particles = parts.n_particles;
    #pragma omp parallel for
    for (int i = 0; i < n_particles; i++){
        acc[i] = calcAcc(i,parts,h,k,n,lambda,nu);
    }
}

#endif