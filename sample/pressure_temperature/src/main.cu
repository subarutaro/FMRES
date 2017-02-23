#include <iostream>
#include <replica_exchanger.hxx>

#include "condition.hxx"
#include "molecular_dynamics.cuh"

static const int ndim = 2;
int main(){
  ReplicaExchanger<ndim,Condition> rem;
  const int size[2] = {1,1};

  rem.Initialize(size);
  double max[ndim] = {2.f,2.f};
  double min[ndim] = {1.0f,1.0f};
  rem.SetConditionRegularInterval(max,min);

  MDParams param;
  param.dt = 0.0005f;
  param.rcut = 4.0f;
  param.nmol = 256;
  param.nreplica = size[0] * size[1];
  param.ngpu = 0;
  param.mass_tstat = 1.f;
  param.mass_bstat = 1.f;

  MultipleMolecularDynamics mmd;
  mmd.Initialize(param);
  mmd.GenerateFCC(0.8f);
  mmd.SetCondition(0,rem.GetCondition(0));
  mmd.GenerateRandomVelocity();

  mmd.isScalingOn = true;
  mmd.isThermostatOn = false;
  mmd.isBarostatOn = true;

  for(int s=-100;s<1000;s++){
    std::cout << s;
    mmd.DebugPrint();
    if(s==0){
      mmd.isScalingOn = false;
      mmd.isThermostatOn = true;
    }
    mmd.MultipleStep(10);
  }
  return 0;
};
