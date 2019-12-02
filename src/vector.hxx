#ifndef HXX_VECTOR
#define HXX_VECTOR

#include <iostream>
#include <iomanip>

template <class T,int D>
class Vector{
 public:
  typedef Vector<T,D> Me;
  #define FORDIM() for(int i=0;i<D;i++)

  T data[D];

  Vector(){};
  Vector(const T& s){FORDIM() data[i] = s;}
  Vector(const Me& v){FORDIM() data[i] = v[i];}
  template <class T2>
  Vector(const Vector<T2,D>& v){FORDIM() data[i] = (T)v[i];}
  ~Vector(){};

  T &operator[](const int index){return data[index];}
  const T &operator[](const int index) const {return data[index];}
  T &operator()(const int index){return data[index];}
  const T &operator()(const int index) const {return data[index];}

  Me operator+(const Me &p) const {
    Me val;
    FORDIM() val[i] = data[i] + p[i];
    return val;
  }
  Me operator-(const Me &p) const {
    Me val;
    FORDIM() val[i] = data[i] - p[i];
    return val;
  }
  Me operator*(const Me &p) const {
    Me val;
    FORDIM() val[i] = data[i] * p[i];
    return val;
  }
  Me operator/(const Me &p) const {
    Me val;
    FORDIM() val[i] = data[i] / p[i];
    return val;
  }
  Me operator+=(const Me &p){
    FORDIM() data[i] += p[i];
    return *this;
  }
  Me operator-=(const Me &p){
    FORDIM() data[i] -= p[i];
    return *this;
  }
  Me operator*=(const Me &p){
    FORDIM() data[i] *= p[i];
    return *this;
  }
  Me operator/=(const Me &p){
    FORDIM() data[i] /= p[i];
    return *this;
  }
  //scalar product
  T operator%(const Me &p) const {
    T val = (T)0;
    FORDIM() val += data[i]*p[i];
    return val;
  }

  Me operator*(Vector<Vector<T,D>,D> &m){
    Me val;
    FORDIM() val[i] = m[i]%(*this);
    return val;
  }

  Me operator+(const T &s) const {
    Me val;
    FORDIM() val[i] = data[i] + s;
    return val;
  }
  Me operator-(const T &s) const {
    Me val;
    FORDIM() val[i] = data[i] - s;
    return val;
  }
  Me operator*(const T &s) const {
    Me val;
    FORDIM() val[i] = data[i] * s;
    return val;
  }
  Me operator/(const T &s) const {
    Me val;
    FORDIM() val[i] = data[i] / s;
    return val;
  }
  Me operator+=(const T &s){
    FORDIM() data[i] += s;
    return *this;
  }
  Me operator-=(const T &s){
    FORDIM() data[i] -= s;
    return *this;
  }
  Me operator*=(const T &s){
    FORDIM() data[i] *= s;
    return *this;
  }
  Me operator/=(const T &s){
    FORDIM() data[i] /= s;
    return *this;
  }

  /*
  template <class T2>
  const Me &operator=(const Vector<T2,D> &p){
    FORDIM() data[i] = p[i];
    return *this;
  }
  //*/
  const Me &operator=(const T &s){
    FORDIM() data[i] = s;
    return *this;
  }

  bool  operator>(const Me &p){
    bool val = true;
    FORDIM() val = val && data[i] >  p[i];
    return val;
  }
  bool  operator>=(const Me &p){
    bool val = true;
    FORDIM() val = val && data[i] >= p[i];
    return val;
  }
  bool  operator<(const Me &p){
    bool val = true;
    FORDIM() val = val && data[i] <  p[i];
    return val;
  }
  bool  operator<=(const Me &p){
    bool val = true;
    FORDIM() val = val && data[i] <= p[i];
    return val;
  }

  //cast operator
  operator       T*()       {return data;}
  operator const T*() const {return data;}
  /*
  template <class P>
  operator Vector<P,D>(){
    Vector<P,D> val;
    FORDIM() val[i] = (P)data[i];
    return val;
  }
  //*/
  // mathmatical operator
  friend Me min(const Me p0,const Me p1){
    Me val;
    FORDIM() val[i] = (p0[i] < p1[i]) ? p0[i]:p1[i];
    return val;
  }
  friend Me max(const Me p0,const Me p1){
    Me val;
    FORDIM() val[i] = (p0[i] > p1[i]) ? p0[i]:p1[i];
    return val;
  }

  // ostream operator
  friend std::ostream& operator<<(std::ostream &os,const Me p){
    FORDIM() os << " " << p[i];
    return os;
  }
#undef FORDIM
};

template <class T,int D>
T sum(Vector<T,D> p){
  T val = 0;
  for(int i=0;i<D;i++) val += p[i];
  return val;
}
template <class T,int D>
T prod(Vector<T,D> p){
  T val = 1;
  for(int i=0;i<D;i++) val *= p[i];
  return val;
}
template <class T,int D>
Vector<T,D> cumprod(Vector<T,D> p){
  Vector<T,D> val;
  T prev = 1;
  for(int i=0;i<D;i++) val[i] = prev = prev*p[i];
  return val;
}
template <class T,int D>
Vector<T,D> inv(Vector<T,D> p){
  Vector<T,D> val;
  for(int i=0;i<D;i++) val[i] = (T)1 / p[i];
  return val;
}

//vector product
template<class T>
Vector<T,3> operator^(const Vector<T,3> &lhs,const Vector<T,3> &rhs){
  Vector<T,3> val;
  val[0] = lhs[1]*rhs[2] - lhs[2]*rhs[1];
  val[1] = lhs[2]*rhs[0] - lhs[0]*rhs[2];
  val[2] = lhs[0]*rhs[1] - lhs[1]*rhs[0];
  return val;
}

template <class T,int N,int M>
Vector<Vector<T,M>,N> Trans(Vector<Vector<T,N>,M> &v){
  Vector<Vector<T,M>,N> tmp;
  for(int i=0;i<N;i++)
    for(int j=0;j<M;j++)
      tmp[i][j] = v[j][i];
  return tmp;
}

template <class T,int N,int M,int O>
Vector<Vector<T,O>,N> operator*(const Vector<Vector<T,M>,N> &lhs,
				const Vector<Vector<T,O>,M> &rhs){
  const Vector<Vector<T,M>,O> rhst = Trans(rhs);
  Vector<Vector<T,O>,N> tmp;
  for(int i=0;i<N;i++)
    for(int j=0;j<O;j++)
      tmp[i][j] = (lhs[i]%rhst[j]);
  return tmp;
}

template <class T,int N,int M>
Vector<T,M> operator*(const Vector<Vector<T,M>,N> &lhs,
		      const Vector<T,M> &rhs){
  Vector<T,N> tmp;
  for(int i=0;i<N;i++) tmp[i] = lhs[i]%rhs;
  return tmp;
}

typedef Vector<double,2> dvec2;
typedef Vector<double,3> dvec3;
typedef Vector<double,4> dvec4;
typedef Vector<int,3>    ivec3;
typedef Vector<int,4>    ivec4;

typedef Vector<dvec3,3> dvec33;
typedef Vector<dvec4,4> dvec44;

#endif
