#ifndef HXX_MOLECULAR_DYNAMICS
#define HXX_MOLECULAR_DYNAMICS

#include <cassert>
//#include <vector_types.h>

#include <algorithm>

#include "condition.hxx"
#include "integrator.cuh"

#include "constant.cuh"

const int MAX_NGPU = 8;

class MDParams{
public:
  union{
    FP data[7];
    struct{
      SP ngpu;
      SP nmol;
      SP nreplica;
      FP dt;
      FP rcut;
      FP mass_tstat;
      FP mass_bstat;
    };
  };
  MDParams(){};
  MDParams(const MDParams& _p){*this = _p;};
  const MDParams& operator=(const MDParams& rhs){
#if 0
    (*this).ngpu        = rhs.ngpu;
    (*this).nmol        = rhs.nmol;
    (*this).nreplica    = rhs.nreplica;
    (*this).dt          = rhs.dt;
    (*this).rcut        = rhs.rcut;
    (*this).mass_tstat  = rhs.mass_tstat;
    (*this).mass_bstat  = rhs.mass_bstat;
#else
    for(int i=0;i<7;i++) data[i] = rhs.data[i];
#endif
    return *this;
  }
};

class MultipleMolecularDynamics{
private:

public:
  MDParams param;
  class Variables{
  private:
  public:
    SP4 *r = nullptr;
    FP4 *v = nullptr;
    FP4 *f = nullptr;

    FP  *s = nullptr;
    FP  *z = nullptr;
    FP4 *L = nullptr; // L.w = volume
    FP  *Ps= nullptr;
    FP  *Pv= nullptr;

    SP4 *d_r[MAX_NGPU];
    FP4 *d_v[MAX_NGPU];
    FP4 *d_f[MAX_NGPU];

    FP  *d_s[MAX_NGPU];
    FP  *d_z[MAX_NGPU];
    FP  *d_L[MAX_NGPU];
    FP  *d_Ps[MAX_NGPU];
    FP  *d_Pv[MAX_NGPU];

    Variables(){};
    ~Variables(){
      if(r != nullptr) delete[] r;
      if(v != nullptr) delete[] v;
      if(f != nullptr) delete[] f;
      if(s != nullptr) delete[] s;
      if(s != nullptr) delete[] z;
      if(L != nullptr) delete[] L;
      if(Ps != nullptr) delete[] Ps;
      if(Pv != nullptr) delete[] Pv;
    }
    void Initialize(const MDParams& param){
      assert(sizeof(FP4) == sizeof(SP4));

      r = new SP4[param.nmol * param.nreplica];
      v = new FP4[param.nmol * param.nreplica];
      f = new FP4[param.nmol * param.nreplica];
      s = new FP[param.nreplica];
      z = new FP[param.nreplica];
      L = new FP4[param.nreplica];
      Ps = new FP[param.nreplica];
      Pv = new FP[param.nreplica];

      for(int rep=0;rep<param.nreplica;rep++){
	z[rep] = 0.f;
	s[rep] = 1.f;
      }

      if(param.ngpu < 1) return;
      // Allocation on device memory
      for(int d=0;d<param.ngpu;d++){
	const int nreplica_gpu = param.nreplica / param.ngpu + (d < (param.nreplica%param.ngpu))? 1 : 0;
	cudaMalloc((void**)&d_r[d], sizeof(FP4)*param.nmol*nreplica_gpu);
	cudaMalloc((void**)&d_v[d], sizeof(FP4)*param.nmol*nreplica_gpu);
	cudaMalloc((void**)&d_f[d], sizeof(FP4)*param.nmol*nreplica_gpu);
	cudaMalloc((void**)&d_s[d], sizeof(FP)*nreplica_gpu);
	cudaMalloc((void**)&d_z[d], sizeof(FP)*nreplica_gpu);
	cudaMalloc((void**)&d_L[d], sizeof(FP)*nreplica_gpu);
	cudaMalloc((void**)&d_Ps[d],sizeof(FP)*nreplica_gpu);
	cudaMalloc((void**)&d_Pv[d],sizeof(FP)*nreplica_gpu);
      }
    }
    void MemcpyH2D(const MDParams& param){
      assert(param.ngpu>0);
      int count = 0;
      for(int d=0;d<param.ngpu;d++){
	const int nreplica_gpu = param.nreplica / param.ngpu + (d < (param.nreplica%param.ngpu))? 1 : 0;
	cudaMemcpy(d_r[d], r+count*param.nmol, sizeof(FP4)*param.nmol*nreplica_gpu,cudaMemcpyHostToDevice);
	cudaMemcpy(d_v[d], v+count*param.nmol, sizeof(FP4)*param.nmol*nreplica_gpu,cudaMemcpyHostToDevice);
	cudaMemcpy(d_f[d], f+count*param.nmol, sizeof(FP4)*param.nmol*nreplica_gpu,cudaMemcpyHostToDevice);
	cudaMemcpy(d_s[d], s+count, sizeof(FP)*nreplica_gpu,cudaMemcpyHostToDevice);
	cudaMemcpy(d_z[d], s+count, sizeof(FP)*nreplica_gpu,cudaMemcpyHostToDevice);
	cudaMemcpy(d_L[d], L+count, sizeof(FP)*nreplica_gpu,cudaMemcpyHostToDevice);
	cudaMemcpy(d_Ps[d],Ps+count,sizeof(FP)*nreplica_gpu,cudaMemcpyHostToDevice);
	cudaMemcpy(d_Pv[d],Pv+count,sizeof(FP)*nreplica_gpu,cudaMemcpyHostToDevice);
	count += nreplica_gpu;
      }
    }
    void MemcpyD2H(const MDParams& param){
      assert(param.ngpu>0);
      int count = 0;
      for(int d=0;d<param.ngpu;d++){
	const int nreplica_gpu = param.nreplica / param.ngpu + (d<(param.nreplica%param.ngpu))? 1 : 0;
	cudaMemcpy(r+count*param.nmol, d_r[d], sizeof(FP4)*param.nmol*nreplica_gpu,cudaMemcpyDeviceToHost);
	cudaMemcpy(v+count*param.nmol, d_v[d], sizeof(FP4)*param.nmol*nreplica_gpu,cudaMemcpyDeviceToHost);
	cudaMemcpy(f+count*param.nmol, d_f[d], sizeof(FP4)*param.nmol*nreplica_gpu,cudaMemcpyDeviceToHost);
	cudaMemcpy(s+count, d_s[d], sizeof(FP)*nreplica_gpu,cudaMemcpyDeviceToHost);
	cudaMemcpy(z+count, d_z[d], sizeof(FP)*nreplica_gpu,cudaMemcpyDeviceToHost);
	cudaMemcpy(L+count, d_L[d], sizeof(FP)*nreplica_gpu,cudaMemcpyDeviceToHost);
	cudaMemcpy(Ps+count,d_Ps[d],sizeof(FP)*nreplica_gpu,cudaMemcpyDeviceToHost);
	cudaMemcpy(Pv+count,d_Pv[d],sizeof(FP)*nreplica_gpu,cudaMemcpyDeviceToHost);
	count += nreplica_gpu;
      }
    }
  };
  Variables vars;

  Condition *cond;
  Result    *result;
  bool isScalingOn   = false;
  bool isThermostatOn = false;
  bool isBarostatOn  = false;

  MultipleMolecularDynamics(){
  }
  ~MultipleMolecularDynamics(){
    if(cond != nullptr) delete[] cond;
    if(result != nullptr) delete[] result;
  }

  void Initialize(const MDParams& _param){
    param = _param;
    vars.Initialize(_param);
    cond = new Condition[param.nreplica];
    result = new Result[param.nreplica];
  }

  template <bool CALC_POT>
  void CalcLJ(){
    for(int rep=0;rep<param.nreplica;rep++){
      const int4 *r = vars.r + rep * param.nmol;
      float4 *f = vars.f + rep * param.nmol;
      const float4 L = vars.L[rep];
      float4 virial = {0.f,0.f,0.f,0.f};
      if(0.5f*std::min({L.x,L.y,L.z}) < param.rcut){
	param.rcut = 0.5f*std::min({L.x,L.y,L.z});
	std::cerr << "warning: rcut > L. rcut is set as 0.5*L(=" << param.rcut << ")" << std::endl;
      }
      const float rcut2 = param.rcut * param.rcut;
      if(CALC_POT) result[rep].vdw = 0.f;
      for(int i=0;i<param.nmol;i++){
	const SP4 r_i = r[i];
	FP4 f_i = {0.f,0.f,0.f,0.f};
	for(int j=0;j<param.nmol;j++){
	  if(i==j) continue;
	  const FP dx = I2F*(FP)(r_i.x - r[j].x)*L.x;
	  const FP dy = I2F*(FP)(r_i.y - r[j].y)*L.y;
	  const FP dz = I2F*(FP)(r_i.z - r[j].z)*L.z;
	  const FP r2 = dx*dx + dy*dy + dz*dz;
	  if(r2 < rcut2){
	    const FP r2i = 1.f / r2;
	    const FP r6i = r2i * r2i * r2i;
	    const FP ftmp = (48.f * r6i - 24.f) * r6i * r2i;
	    f_i.x += ftmp * dx;
	    f_i.y += ftmp * dy;
	    f_i.z += ftmp * dz;
	    virial.x += ftmp*dx*dx;
	    virial.y += ftmp*dy*dy;
	    virial.z += ftmp*dz*dz;
	    if(CALC_POT) f_i.w += 4.f * r6i * (r6i - 1.f);
	  }
	}
	f[i] = f_i;
	if(CALC_POT) result[rep].vdw += 0.5f*f_i.w;
      }
      virial.x *= 0.5f;
      virial.y *= 0.5f;
      virial.z *= 0.5f;
      virial.w = virial.x + virial.y + virial.z;
      result[rep].virial = virial;
    }
  }
  void CalcKin(){
    for(int rep=0;rep<param.nreplica;rep++){
      const FP4 *v = vars.v + rep*param.nmol;
      const FP4 L = vars.L[rep];
      float kin = 0.f;
      for(int i=0;i<param.nmol;i++){
	kin += v[i].w * v[i].x*v[i].x/(L.x*L.x);
	kin += v[i].w * v[i].y*v[i].y/(L.y*L.y);
	kin += v[i].w * v[i].z*v[i].z/(L.z*L.z);
      }
      result[rep].kinetic = 0.5f*kin;
    }
  }
  void CalcHamiltonian(){
    for(int rep=0;rep<param.nreplica;rep++){
      Result &r = result[rep];
      r.hamiltonian = 0.f;
      r.hamiltonian += r.vdw + r.kinetic;
      if(isThermostatOn){
	r.hamiltonian += 3.f * param.nmol * cond[rep].temperature * log(vars.s[rep]);
	r.hamiltonian += 0.5f * vars.z[rep] * vars.z[rep] * param.mass_tstat;
      }
      if(isBarostatOn){
	r.hamiltonian += cond[rep].pressure * vars.L[rep].w;
	r.hamiltonian += 0.5f * vars.Pv[rep] * vars.Pv[rep] * param.mass_bstat;
      }
    }
  }

  void GenerateFCC(const float density){
    //assert(isInitialized);
    int n=1;
    while(4*n*n*n<param.nmol) n++;
    const FP unit = 1.f / (FP)n;
    const FP dx[4] = {0.f, 0.5f*unit, 0.f,       0.5f*unit};
    const FP dy[4] = {0.f, 0.f,       0.5f*unit, 0.5f*unit};
    const FP dz[4] = {0.f, 0.5f*unit, 0.5f*unit, 0.f};
    for(int rep=0;rep<param.nreplica;rep++){
      int count = 0;
      SP4 *r = vars.r + rep*param.nmol;
      for(int z=0;z<n;z++){
	for(int y=0;y<n;y++){
	  for(int x=0;x<n;x++){
	    for(int i=0;i<4;i++){
	      if(count==param.nmol) break;
	      r[count].x = (SP)(F2I*(unit*x + dx[i] - 0.5f));
	      r[count].y = (SP)(F2I*(unit*y + dy[i] - 0.5f));
	      r[count].z = (SP)(F2I*(unit*z + dz[i] - 0.5f));
	      r[count].w = count;
	      //std::cout << r[count].x << ' ' << r[count].y << ' ' << r[count].y << std::endl;
	      count++;
	    }
	  }
	}
      }
      assert(count==param.nmol);
      FP4 &L = vars.L[rep];
      L.x =  L.y = L.z = powf((FP)param.nmol/density,1.f/3.f);
      L.w = L.x * L.y * L.z;
    }
  }
  void GenerateRandomVelocity(){
    for(int rep=0;rep<param.nreplica;rep++){
      FP4 *v = vars.v + rep*param.nmol;
      for(int i=0;i<param.nmol;i++){
	v[i].x = (FP)rand() / ((FP)RAND_MAX + 1) - 0.5f;
	v[i].y = (FP)rand() / ((FP)RAND_MAX + 1) - 0.5f;
	v[i].z = (FP)rand() / ((FP)RAND_MAX + 1) - 0.5f;
	v[i].w = 1.f; // set mass of atom as 1.0
      }
      RemoveTotalMomentum();
      CalcKin();
      VelocityScaling();
    }
  }

  void RemoveTotalMomentum(){
    for(int rep=0;rep<param.nreplica;rep++){
      FP4 *v = vars.v + rep*param.nmol;
      FP4 sum = {0.f,0.f,0.f,0.f};
      for(int i=0;i<param.nmol;i++){
    	sum.x += v[i].w * v[i].x;
	sum.y += v[i].w * v[i].y;
	sum.z += v[i].w * v[i].z;
	sum.w += v[i].w;
      }
      sum.x /= sum.w;
      sum.y /= sum.w;
      sum.z /= sum.w;
      for(int i=0;i<param.nmol;i++){
	v[i].x -= sum.x;
	v[i].y -= sum.y;
	v[i].z -= sum.z;
      }
    }
  }

  void VelocityScaling(){
    CalcKin();
    for(int rep=0;rep<param.nreplica;rep++){
      const FP scale = sqrt(1.5f*param.nmol*cond[rep].temperature / result[rep].kinetic);
      FP4 *v = vars.v + rep*param.nmol;
      for(int i=0;i<param.nmol;i++){
	v[i].x *= scale;
	v[i].y *= scale;
	v[i].z *= scale;
      }
    }
    RemoveTotalMomentum();
  }

  void OneStep(){
    for(int rep=0;rep<param.nreplica;rep++){
      const int offset = rep * param.nmol;
      SP4 *r = vars.r + offset;
      FP4 *v = vars.v + offset;
      FP4 *f = vars.f + offset;
      FP4 &L = vars.L[rep];
      FP  &z = vars.z[rep];
      FP  &s = vars.s[rep];
      FP  &Pv= vars.Pv[rep];
      const float4 vir = result[rep].virial;
      const float P = cond[rep].pressure;
      if(isBarostatOn) andersen_vel_vir(Pv,vir.w,s,L.w,P,param.dt);
      for(int i=0;i<param.nmol;i++){
	if(isThermostatOn) nose_hoover_vel(v[i],z,param.dt);
	vel_velret_vel(v[i],f[i],L,param.dt);
	vel_velret_pos(r[i],v[i],L,param.dt);
      }
      andersen_volume(L,Pv,param.mass_bstat,s,param.dt);
    }
    CalcLJ<false>();
    CalcKin();
    for(int rep=0;rep<param.nreplica;rep++){
      const int offset = rep * param.nmol;
      FP4 *v = vars.v + offset;
      FP4 *f = vars.f + offset;
      FP4 &L = vars.L[rep];
      FP  &z = vars.z[rep];
      FP  &s = vars.s[rep];
      FP  &Pv= vars.Pv[rep];
      const float kin = result[rep].kinetic;
      const float4 vir = result[rep].virial;
      const float T = cond[rep].temperature;
      const float P = cond[rep].pressure;

      if(isThermostatOn) nose_hoover_zeta(z,s,kin,T,param.mass_tstat,param.dt,param.nmol);
      if(isBarostatOn) andersen_vel_kin(Pv,kin,s,L.w,param.dt);
      for(int i=0;i<param.nmol;i++){
	vel_velret_vel(v[i],f[i],L,param.dt);
	if(isThermostatOn) nose_hoover_vel(v[i],z,param.dt);
      }
      if(isBarostatOn) andersen_vel_vir(Pv,vir.w,s,L.w,P,param.dt);
    }
  }
  void CalcPressure(){
    for(int rep=0;rep<param.nreplica;rep++){
      result[rep].pressure = (result[rep].virial.w + 2.f*result[rep].kinetic)/(3.f*vars.L[rep].w);
    }
  }
  void MultipleStep(const int nstep){
    for(int s=0;s<nstep;s++){
      OneStep();
      if(isScalingOn) VelocityScaling();
    }
    CalcLJ<true>();
    CalcKin();
    CalcPressure();
    CalcHamiltonian();
  }

  void SetCondition(const int i,const Condition& c){cond[i] = c;}
  Result* GetResultPointer(){return result;}

  void DebugPrint(){
    for(int rep=0;rep<param.nreplica;rep++){
#if 1
      const Result r = result[rep];
      std::cout << ' ' << rep;
      std::cout << ' ' << r.vdw;
      std::cout << ' ' << r.kinetic;
      std::cout << ' ' << vars.s[rep];
      std::cout << ' ' << vars.z[rep];
      std::cout << ' ' << vars.L[rep].w;
      std::cout << ' ' << vars.Pv[rep];
      std::cout << ' ' << r.pressure;
      std::cout << ' ' << r.hamiltonian;
      std::cout << std::endl;
#else
      const SP4* r = vars.r + rep*param.nmol;
      const FP4* v = vars.v + rep*param.nmol;
      const FP4 L = vars.L[rep];
      for(int i=0;i<param.nmol;i++){
	const FP4 rf = {I2F*(FP)r[i].x*L.x, I2F*(FP)r[i].y*L.y, I2F*(FP)r[i].z*L.z};
	std::cout << ' ' << i;
	std::cout << ' ' << rf.x   << ' ' << rf.y   << ' ' << rf.z;
	std::cout << ' ' << v[i].x << ' ' << v[i].y << ' ' << v[i].z;
	std::cout << std::endl;
      }
#endif
    }
  }
};

#endif
