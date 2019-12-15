#ifndef HXX_CONDITION
#define HXX_CONDITION

#include <iostream>
#include <sum_log.hxx>
#include <cmath>

#include <vector.hxx>

typedef struct{
  float x;
  float y;
}float2;
typedef struct{
  float x;
  float y;
  float z;
}float3;
typedef struct{
  float x;
  float y;
  float z;
  float w;
}float4;

int NMOL = 108;

class Condition{
public:
  union{
    float data[2];
    struct{
      float temperature;
      float pressure;
    };
  };
  int nmol;

  Condition(){ nmol = NMOL;}

  friend std::ostream& operator<<(std::ostream &os,const Condition& c){
    os << c.temperature << " " << c.pressure;
    return os;
  }
  friend std::istream& operator>>(std::istream &is,Condition &c){
    is >> c.nmol >> c.temperature >> c.pressure;
    NMOL = c.nmol;
    return is;
  }

  float& operator[](const int i){return data[i];}
  const float& operator[](const int i) const {return data[i];}

  float calcLogBoltzmannFactor(const Vector<double,2>& bin) const {
    return -(bin[0] + pressure*nmol/bin[1])/temperature;
  }
};

#define NDATA_RESULT 4
class Result{
public:
  union{
    float data[NDATA_RESULT];
    struct{
      float pot;
      float order_parameter;
      float dens;
      float smectic_parameter;
    };
  };
  Result(){for(int i=0;i<NDATA_RESULT;i++)data[i] = 0.f;}
  ~Result(){}

  Vector<double,2> getBin() const {
    Vector<double,2> ret;
    ret[0] = pot;
    ret[1] = dens;
    return ret;
  }

  float& operator[](const int i){ return data[i];}
  const float& operator[](const int i) const { return data[i];}

  Result operator+(const Result& rhs) const {
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++) tmp[i] = data[i] + rhs[i];
    return tmp;
  }
  Result operator+=(const Result& rhs){
    *this = *this + rhs;
    return *this;
  }
  Result operator-(const Result& rhs) const {
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++) tmp[i] = data[i] - rhs[i];
    return tmp;
  }
  Result operator-=(const Result& rhs){
    *this = *this - rhs;
    return *this;
  }

  Result operator*(const Result& rhs) const {
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++) tmp[i] = data[i] * rhs[i];
    return tmp;
  }
  Result operator=(const Result& rhs){
    for(int i=0;i<NDATA_RESULT;i++) data[i] = rhs[i];
    return *this;
  }

  Result operator+(const float& rhs){
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++) tmp[i] = data[i] + rhs;
    return tmp;
  }
  Result operator-(const float& rhs){
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++) tmp[i] = data[i] - rhs;
    return tmp;
  }
  Result operator/(const double& rhs){
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++) tmp[i] = data[i] / rhs;
    return tmp;
  }
  Result operator=(const float& rhs){
    for(int i=0;i<NDATA_RESULT;i++) data[i] = rhs;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream &os,const Result &r){
    for(int i=0;i<NDATA_RESULT;i++) os << " " << r.data[i];
    return os;
  }
  friend std::istream& operator>>(std::istream &is,Result &r){
    int step;
    is >> step;
    //>> r.pot >> r.order_parameter >> r.dens >> r.smectic_parameter;
    for(int i=0;i<NDATA_RESULT;i++) is >> r.data[i];
    return is;
  }
  friend Result log(const Result& r){
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++) tmp[i] = log(r[i]);
    return tmp;
  }
  friend Result exp(const Result& r){
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++) tmp[i] = exp(r[i]);
    return tmp;
  }
  friend Result fabs(const Result& r){
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++) tmp[i] = fabs(r[i]);
    return tmp;
  }
  friend Result max(const Result& a,const Result& b){
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++)
      tmp[i] = (a[i]>b[i]) ? a[i] : b[i];
    return tmp;
  }
  friend Result min(const Result& a,const Result& b){
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++)
      tmp[i] = (a[i]<b[i]) ? a[i] : b[i];
    return tmp;
  }
  friend Result sum_log(const Result& a,const Result& b){
    Result tmp;
    for(int i=0;i<NDATA_RESULT;i++) tmp[i] = sum_log(a[i],b[i]);
    return tmp;
  }
};


#endif // HXX_CONDITION
