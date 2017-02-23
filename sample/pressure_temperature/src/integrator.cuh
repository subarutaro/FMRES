#ifndef CUH_INTEGRATOR
#define CUH_INTEGRATOR

#include "constant.cuh"

__forceinline__
__device__ __host__
void vel_velret_pos(int4 &r,const float4 v,const float4 cell_size,const float dt){
  r.x += (SP)(F2I * v.x / (cell_size.x * cell_size.x)* dt);
  r.y += (SP)(F2I * v.y / (cell_size.y * cell_size.y)* dt);
  r.z += (SP)(F2I * v.z / (cell_size.z * cell_size.z)* dt);
}

__forceinline__
__device__ __host__
void vel_velret_vel(float4 &v,const float4 force,const float4 cell_size,const float dt){
  const float massi = 1.f / v.w;
  v.x += force.x * massi * cell_size.x * 0.5f * dt;
  v.y += force.y * massi * cell_size.y * 0.5f * dt;
  v.z += force.z * massi * cell_size.z * 0.5f * dt;
}

__forceinline__
__device__ __host__
void nose_hoover_vel(float4 &v, const float z,const float dt){
  const float tmp = exp(-0.5f*z*dt);
  v.x *= tmp;
  v.y *= tmp;
  v.z *= tmp;
}

__forceinline__
__device__ __host__
void nose_hoover_zeta
(float &z, float &s,
 const float kin,const float temp,const float Q,const float dt,const int nmol){
  z += (2.f*kin - 3.f*nmol*temp)/Q*dt;
  s *= exp(z*dt);
}

__forceinline__
__device__ __host__
void andersen_vel_kin(float &Pv,const float kin,const float s,const float V,const float dt){
  Pv += 2.f * kin / (3.f*s*V) * dt;
}

__forceinline__
__device__ __host__
void andersen_vel_vir(float &Pv,const float vir,const float s,const float V,const float P,const float dt){
  Pv += s*(vir/(3.f*V) - P) * 0.5f*dt;
}

__forceinline__
__device__ __host__
void andersen_volume(float4 &L,const float Pv,const float M,const float s,const float dt){
  L.w += s*Pv/M*0.5f*dt;
  L.x = L.y = L.z = powf(L.w,1.f/3.f);
}
#endif // CUH_INTEGRATOR
