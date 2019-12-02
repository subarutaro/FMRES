#ifndef HXX_SUM_LOG
#define HXX_SUM_LOG

#include <cmath>

template <class T>
T max(const T a,const T b){
  return (a>b)?a:b;
}
template <class T>
T min(const T a,const T b){
  return (a<b)?a:b;
}

template <class T>
inline T sum_log(const T a,const T b){
  if(a==0.0 || b==0.0) return a+b;
  return max(a,b) + log(exp(-fabs(b-a))+1.);
}

#endif
