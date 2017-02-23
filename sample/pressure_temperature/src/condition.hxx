#ifndef HXX_CONDITION
#define HXX_CONDITION

#include <iostream>
#include <sstream>
#include <iomanip>

#include <string>

#include <cassert>

class Condition{
public:
  union{
    float data[2];
    struct{
      float pressure;
      float temperature;
    };
  };

  float& operator [](const int i){
    assert(i<2);
    return data[i];
  }
  const float& operator [](const int i) const {
    assert(i<2);
    return data[i];
  }

  const Condition operator=(Condition rhs){
    (*this).pressure = rhs.pressure;
    (*this).temperature = rhs.temperature;
    return *this;
  }

  const std::string FileName(){
    std::stringstream strs;
    strs << 'P';
    strs << std::fixed << std::setprecision(8);
    strs << pressure;
    strs << 'T';
    strs << std::fixed << std::setprecision(8);
    strs << temperature;
    strs << ".dat";
    return strs.str();
  }
};

class Result{
public:
  union{
    float data[14];
    struct{
      float vdw;
      float coulomb;
      float wall;
      float kinetic;
      float volume;
      float4 virial;
      float pressure;
      float temperature;
      float eng_tstat;
      float eng_bstat;
      float hamiltonian;
    };
  };

  const Result& operator=(const Result& rhs){
    for(int i=0;i<14;i++) data[i] = rhs.data[i];
    return (*this);
  }
};

#endif
