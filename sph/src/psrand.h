#ifndef _PSRAND_H_
#define _PSRAND_H_

#include <stdlib.h>
#include "log.h"

/*sets random seed*/
static inline void set_seed(int seed){
    LOG("set seed = %d",seed);
    srand(seed);
}

/*returns a random float in `[0,1]`*/
static inline float frand(void){
    return (float)rand()/(float)RAND_MAX;
}

/*returns a random float in `[rmin,rmax]`*/
static inline float frand_range(float rmin, float rmax){
    float r = frand();
    float range = rmax - rmin;
    return rmin + range * r;
}

#endif