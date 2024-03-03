#ifndef _VECTORS_H_
#define _VECTORS_H_

#include <math.h>

typedef struct vec3{
    float x,y,z;
} vec3;

static inline vec3 v3(float x, float y, float z){
    vec3 out;
    out.x = x;
    out.y = y;
    out.z = z;
    return out;
}

static inline vec3 v3v3add(vec3 l, vec3 r){
    vec3 out;
    out.x = l.x + r.x;
    out.y = l.y + r.y;
    out.z = l.z + r.z;
    return out;
}

static inline vec3 v3v3sub(vec3 l, vec3 r){
    vec3 out;
    out.x = l.x - r.x;
    out.y = l.y - r.y;
    out.z = l.z - r.z;
    return out;
}

static inline vec3 v3v3mul(vec3 l, vec3 r){
    vec3 out;
    out.x = l.x * r.x;
    out.y = l.y * r.y;
    out.z = l.z * r.z;
    return out;
}

static inline vec3 v3v3div(vec3 l, vec3 r){
    vec3 out;
    out.x = l.x / r.x;
    out.y = l.y / r.y;
    out.z = l.z / r.z;
    return out;
}

static inline float v3len2(vec3 l){
    return l.x*l.x + l.y*l.y + l.z*l.z;
}

static inline float v3len(vec3 l){
    return sqrtf(v3len2(l));
}

static inline vec3 v3fmul(vec3 l, float r){
    vec3 out;
    out.x = l.x * r;
    out.y = l.y * r;
    out.z = l.z * r;
    return out;
}

static inline vec3 v3zero(void){
    vec3 out;
    out.x = 0;
    out.y = 0;
    out.z = 0;
    return out;
}

/*static inline Vector3 vec3toVector3(vec3 v){
    Vector3 out;
    out.x = v.x;
    out.y = v.y;
    out.z = v.z;
    return out;
}*/

#endif