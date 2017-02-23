#include <cstdio>
#include <cstdlib>
#include <cassert>

template<class ...Args>
int calc_size(Args args...,const int i,const int size,const int sum){
  return calc_size(sum+i*size,args);
}

template<>
int calc_size(const int sum){return sum;}

template <int N,class T>
class NDArray{
private:
  T* data = nullptr;
  int size[N];
public:
  void allocate(const int _size[N]){
    int s = 1;
    for(int i=0;i<N;i++){
      assert(_size[i]>0);
      s *= _size[i];
      size[i] = _size[i];
    }
    if(data != nullptr) free(data);
    if((data = (T*)malloc(s*sizeof(T))) == NULL){
      fprintf(stderr,"error: malloc data in NDimensionalArray failed\n");
      exit(EXIT_FAILURE);
    }
  }
  template <int M,class ...Args>
  T& at(int i, Args args...){data[calc_size(args,i,size[M-1],0)];}
};
