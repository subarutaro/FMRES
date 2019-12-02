#ifndef HXX_ARRAY
#define HXX_ARRAY

#include <vector.hxx>

template<class T,int D>
class Array{
public:
  T* data = nullptr;
  Vector<int,D> sizes;
  Vector<int,D> offsets;
  int Index(Vector<int,D> index) const { return index%offsets;}

  Array():data(nullptr),sizes(0),offsets(Vector<int,D>(0)){}
  Array(const Vector<int,D> &_sizes):data(new T[prod(_sizes)]),sizes(_sizes),offsets(cumprod(_sizes)/_sizes){}
  Array(const Array<T,D> &a){ *this = a;}

  ~Array(){if(data!=nullptr) delete[] data;}

  void Allocate(Vector<int,D> _sizes){
    sizes = _sizes;
    offsets = cumprod(sizes) / sizes;
    if(data!=nullptr) delete[] data;
    data = new T[prod(_sizes)];
  }

  const Vector<int,D> Sizes() const {return sizes;}
  const Vector<int,D> Offsets() const {return offsets;}
  const T& At(const Vector<T,D> index) const { return data[Index(index)];}
  const T& operator()(const Vector<int,D> index) const {return At(index);} // point access,right value
  T& operator()(const Vector<int,D> index){return data[Index(index)];}     // point access, left value

  const T& operator[](const int index) const {return data[index];} // lenear access
  T& operator[](const int index){return data[index];} // lenear access

  operator T*(){return data;}
  operator const T*() const {return data;}
};

#endif
