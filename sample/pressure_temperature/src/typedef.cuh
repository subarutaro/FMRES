#ifndef CUH_TYPEDEF
#define CUH_TYPEDEF

//typedef FP double;
typedef float  FP;
typedef float3 FP3;
typedef float4 FP4;
typedef int    SP;
typedef int3   SP3;
typedef int4   SP4;


__device__ __host__
FP3 operator+(const FP3& lhs,const FP3& rhs){
  return make_float3(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
}
__device__ __host__
FP3 operator-(const FP3& lhs,const FP3& rhs){
  return make_float3(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
}
__device__ __host__
FP3 operator*(const FP3& lhs,const FP3& rhs){
  return make_float3(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
}
__device__ __host__
FP3 operator*(const FP3& lhs,const FP& rhs){
  return make_float3(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
}
__device__ __host__
FP3 operator*(const FP& lhs,const FP3& rhs){
  return rhs*lhs;
}
__device__ __host__
FP3 operator/(const FP3& lhs,const FP3& rhs){
  return make_float3(lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z);
}
__device__ __host__
FP3 operator/(const FP3& lhs,const FP& rhs){
  return make_float3(lhs.x/rhs, lhs.y/rhs, lhs.z/rhs);
}

__device__ __host__
FP3 operator*=(FP3& lhs,const FP& rhs){
  lhs = lhs*rhs;
  return lhs;
}
__device__ __host__
FP3 operator/=(FP3& lhs,const FP& rhs){
  lhs = lhs/rhs;
  return lhs;
}
//inner product
__device__ __host__
FP operator%(const FP3& lhs,const FP3& rhs){
  return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
}
//outer product
__device__ __host__
FP3 operator^(const FP3& lhs,const FP3& rhs){
  return make_float3(lhs.y*rhs.z - lhs.z*rhs.y,
		     lhs.z*rhs.x - lhs.x*rhs.z,
		     lhs.x*rhs.y - lhs.y*rhs.x);
}
__device__ __host__
FP norm2(FP3 v){return v%v;};
__device__ __host__
FP norm(FP3 v){return sqrtf(v%v);};

__device__ __host__
FP4 operator+(const FP4& lhs,const FP4& rhs){
  return make_float4(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z, lhs.w+rhs.w);
}
__device__ __host__
FP4 operator-(const FP4& lhs,const FP4& rhs){
  return make_float4(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z, lhs.w-rhs.w);
}
__device__ __host__
FP4 operator*(const FP4& lhs,const FP4& rhs){
  return make_float4(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z, lhs.w*rhs.w);
}
__device__ __host__
FP4 operator*(const FP4& lhs,const FP& rhs){
  return make_float4(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs, lhs.w*rhs);
}
__device__ __host__
FP4 operator*(const FP& lhs,const FP4& rhs){
  return rhs*lhs;
}
__device__ __host__
FP4 operator/(const FP4& lhs,const FP4& rhs){
  return make_float4(lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z, lhs.w/rhs.w);
}
__device__ __host__
FP4 operator/(const FP4& lhs,const FP& rhs){
  return make_float4(lhs.x/rhs, lhs.y/rhs, lhs.z/rhs, lhs.w/rhs);
}
//inner product
__device__ __host__
FP operator%(const FP4& lhs,const FP4& rhs){
  return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z + lhs.w*rhs.w;
}

class Tensor{
public:
  union{
    FP data[9];
    struct{
      FP xx,xy,xz;
      FP yx,yy,yz;
      FP zx,zy,zz;
    };
  };

  __device__ __host__
  Tensor(){}
  __device__ __host__
  Tensor(const Tensor& _tensor){*this = _tensor;}
  __device__ __host__
  ~Tensor(){}

  __device__ __host__
  const Tensor& operator=(const Tensor& rhs){
    for(int i=0;i<9;i++) data[i] = rhs[i];
    return (*this);
  }
  __device__ __host__
  FP& operator[](const int i){return data[i];};
  __device__ __host__
  const FP& operator[](const int i) const {return data[i];};

  __device__ __host__
  const Tensor operator+(const Tensor& rhs){
    Tensor tmp;
    for(int i=0;i<9;i++) tmp[i] = data[i] + rhs[i];
    return tmp;
  }
  /*
  const Tensor& operator+(const FP& rhs){
    Tensor tmp;
    for(int i=0;i<9;i++) tmp[i] = data[i] + rhs;
    return tmp;
  }
  //*/
  __device__ __host__
  const Tensor operator-(const Tensor& rhs){
    Tensor tmp;
    for(int i=0;i<9;i++) tmp[i] = data[i] - rhs[i];
    return tmp;
  }
  /*
  const Tensor& operator-(const FP& rhs){
    Tensor tmp;
    for(int i=0;i<9;i++) tmp[i] = lhs[i] - rhs;
    return tmp;
  }
  //*/

  //scalar product
  __device__ __host__
  const Tensor operator*(const FP& rhs){
    Tensor tmp;
    for(int i=0;i<9;i++) tmp[i] = data[i] * rhs;
    return tmp;
  }
  //inner product
  __device__ __host__
  const Tensor operator*(const Tensor& rhs){
    Tensor tmp;
    tmp.xx = xx*rhs.xx + xy*rhs.yx + xz*rhs.zx;
    tmp.xy = xx*rhs.xy + xy*rhs.yy + xz*rhs.zy;
    tmp.xz = xx*rhs.xz + xy*rhs.yz + xz*rhs.zz;

    tmp.yx = yx*rhs.xx + yy*rhs.yx + yz*rhs.zx;
    tmp.yy = yx*rhs.xy + yy*rhs.yy + yz*rhs.zy;
    tmp.yz = yx*rhs.xz + yy*rhs.yz + yz*rhs.zz;

    tmp.zx = zx*rhs.xx + zy*rhs.yx + zz*rhs.zx;
    tmp.zy = zx*rhs.xy + zy*rhs.yy + zz*rhs.zy;
    tmp.zz = zx*rhs.xz + zy*rhs.yz + zz*rhs.zz;
    return tmp;
  }

  __device__ __host__
  friend const FP3 operator*(const FP3& lhs,const Tensor& rhs){
    FP3 tmp;
    tmp.x = lhs.x*rhs.xx + lhs.y*rhs.yx + lhs.z*rhs.zx;
    tmp.y = lhs.x*rhs.xy + lhs.y*rhs.yy + lhs.z*rhs.zy;
    tmp.z = lhs.x*rhs.xz + lhs.y*rhs.yz + lhs.z*rhs.zz;
    return tmp;
  }
  __device__ __host__
  friend const FP3 operator*(const FP4& lhs,const Tensor& rhs){
    FP3 tmp = {lhs.x,lhs.y,lhs.z};
    tmp = tmp * rhs;
    return tmp;
  }
  __device__ __host__
  friend const FP3 operator*(const Tensor& lhs,const FP3& rhs){
    FP3 tmp;
    tmp.x = lhs.xx*rhs.x + lhs.xy*rhs.y + lhs.xz*rhs.z;
    tmp.y = lhs.yx*rhs.x + lhs.yy*rhs.y + lhs.yz*rhs.z;
    tmp.z = lhs.zx*rhs.x + lhs.zy*rhs.y + lhs.zz*rhs.z;
    return tmp;
  }
  __device__ __host__
  friend const FP3 operator*(const Tensor& lhs,const FP4& rhs){
    FP3 tmp = {rhs.x,rhs.y,rhs.z};
    tmp = lhs*tmp;
    return tmp;
  }

  __device__ __host__
  const Tensor operator/(const FP& rhs){
    Tensor tmp;
    FP rhsi = 1.f / rhs;
    return (*this) * rhsi;
  }

  /*
  //tensor product
  const Tensor& operator^(const Tensor& rhs){
    Tensor tmp;
    tmp.xx = 

    return tmp;
  }
  //*/
  __device__ __host__
  const Tensor& operator+=(const Tensor& rhs){
    for(int i=0;i<9;i++) data[i] += rhs[i];
    return *this;
  };
  __device__ __host__
  const Tensor& operator-=(const Tensor& rhs){
    for(int i=0;i<9;i++) data[i] -= rhs[i];
    return *this;
  };
};

__device__ __host__
Tensor Trans(const Tensor& t){
  Tensor tmp;
  tmp.xx = t.xx; tmp.xy = t.yx; tmp.xz = t.zx;
  tmp.yx = t.xy; tmp.yy = t.yy; tmp.yz = t.zy;
  tmp.zx = t.xz; tmp.zy = t.yz; tmp.zz = t.zz;
  return tmp;
}

__device__ __host__
FP Trace(const Tensor& t){
  return t.xx + t.yy + t.zz;
}
__device__ __host__
FP Det(const Tensor& t){
  FP tmp = 0.f;
  tmp += t.xx * t.yy * t.zz;
  tmp += t.yx * t.zy * t.xz;
  tmp += t.zx * t.xy * t.yz;
  tmp -= t.xx * t.zy * t.zx;
  tmp -= t.zx * t.yy * t.xz;
  tmp -= t.yx * t.xy * t.zz;
  return tmp;
}

__device__ __host__
Tensor Inverse(const Tensor& t){
  Tensor tmp;
  tmp.xx = t.yy*t.zz - t.yz*t.zy;
  tmp.xy = t.xz*t.zy - t.xy*t.zz;
  tmp.xz = t.xy*t.yz - t.xz*t.yy;

  tmp.yx = t.yz*t.zx - t.xy*t.zz;
  tmp.yy = t.xx*t.zz - t.xz*t.zx;
  tmp.yz = t.xz*t.yx - t.xx*t.yz;

  tmp.zx = t.yx*t.zy - t.yy*t.zx;
  tmp.zy = t.xy*t.zx - t.xx*t.zy;
  tmp.zz = t.xx*t.yy - t.xy*t.yx;
  return tmp / Det(t);
}

#endif
