#include <cassert>
int factorial(const int n){
  assert(n>=0);
  if(n==0) return 1;
  return n*factorial(n-1);
}

int combination(const int a,const int b){
  return factorial(a) / (factorial(a-b) * factorial(b));
}

double random_number(){
  return (double)rand() / ((double)RAND_MAX + 1.0);
}

class Permutation{
private:
  const int  size;
  int  n;
  int  *digit  = nullptr;
  bool *isUsed = nullptr;

  int  **list = nullptr;

  void Permute(int p){
    if(p == size){
      for(int i=0;i<size;i++){
	list[n][i] = digit[i];
      }
      n++;
      return;
    }
    for(digit[p]=0;digit[p]<size;digit[p]++){
      if(isUsed[digit[p]]) continue;
      isUsed[digit[p]] = true;
      Permute(p+1);
      isUsed[digit[p]] = false;
    }
  }
public:
  Permutation(const int _size):size(_size){
    if(digit  == nullptr) digit  = new int[size];
    if(isUsed == nullptr) isUsed = new bool[size];
    if(list == nullptr){
      list = new int*[factorial(size)];
      for(int i=0;i<factorial(size);i++) list[i] = new int[size];
    }
    n = 0;
    for(int i=0;i<size;i++){
      isUsed[i] = false;
    }
    Permute(0);
  }
  ~Permutation(){
    if(digit!=nullptr)  delete[] digit;
    if(isUsed!=nullptr) delete[] isUsed;
    if(list!=nullptr){
      for(int i=0;i<size;i++) delete[] list[i];
      delete[] list;
    }
  }
  const int* operator[](const int i) const {
    assert(i<factorial(size));
    return list[i];
  }
};
