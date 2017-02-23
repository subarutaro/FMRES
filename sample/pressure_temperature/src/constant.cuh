#ifndef CUH_CONSTANT
#define CUH_CONSTANT

//typedef FP double;
typedef float  FP;
typedef float3 FP3;
typedef float4 FP4;
typedef int    SP;
typedef int3   SP3;
typedef int4   SP4;

const FP I2F = 1.f/(FP)(1LL<<32);
const FP F2I = (FP)(1LL<<32);

#endif
