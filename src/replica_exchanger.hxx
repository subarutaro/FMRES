#include <cmath>
#include <cassert>

#include <iostream>
#include <fstream>

#include <string>

#include "permutation.hxx"

template <class T>
void swap(T& a,T& b){
  T tmp = b;
  b = a;
  a = tmp;
}

class AcceptRatio{
private:
public:
  long accept = 0;
  long refuse = 0;
  void Accept(){accept++;}
  void Refuse(){refuse++;}
  AcceptRatio(){};
  ~AcceptRatio(){};
  double Get(){
    if((accept+refuse)==0) return -1;
    return (double)accept / (double)(accept+refuse);
  }
  void Flush(){accept = 0; refuse = 0;};
};

template <int NDIM,class TCondition>
class ReplicaExchanger{
private:
  int size[NDIM];
  int interval[NDIM];
  int nreplica;
  TCondition *cond = nullptr;
  int *index = nullptr;
  bool *isExchanged = nullptr;
  bool isInitialized = false;

  AcceptRatio *accept_ratio = nullptr;

  int dim;
  bool isEven;
  
  std::fstream *output = nullptr;
public:
  ReplicaExchanger(){
    for(int i=0;i<NDIM;i++) size[i] = 1;
    nreplica = 1;
  };
  ~ReplicaExchanger(){
    if(cond  != nullptr) delete[] cond;
    if(index != nullptr) delete[] index;
  }

  void Initialize(const int _size[NDIM]){
    assert(!isInitialized);
    for(int i=0;i<NDIM;i++){
      assert(_size[i]>0);
      size[i] = _size[i];
    }
    nreplica = 1;
    for(int i=0;i<NDIM;i++) nreplica *= size[i];

    if(cond==nullptr)  cond  = new TCondition[nreplica];
    if(index==nullptr) index = new int[nreplica];
    for(int i=0;i<nreplica;i++) index[i] = i;
    if(isExchanged==nullptr) isExchanged = new bool[nreplica];
    for(int i=0;i<nreplica;i++) isExchanged[i] = false;
    if(accept_ratio==nullptr) accept_ratio = new AcceptRatio[NDIM*nreplica];

    for(int d=0;d<NDIM;d++){
      interval[d] = 1;
      for(int dd=0;dd<d;dd++){
	interval[d] *= size[dd];
      }
    }

    dim = 0;
    isEven = true;

    isInitialized = true;
  };

  void SetConditionRegularInterval(const double max[NDIM],const double min[NDIM]){
    double diff[NDIM];
    for(int i=0;i<NDIM;i++) diff[i] = (size[i]>1)?(max[i] - min[i]) / (double)(size[i]-1) : min[i];
    for(int i=0;i<nreplica;i++){
      int ind[NDIM];
      for(int d=0;d<NDIM;d++){
	ind[d] = (i/interval[d])%size[d];
      }
      for(int d=0;d<NDIM;d++) cond[i][d] = min[d] + diff[d]*ind[d];
    }
  }
  /*
  void SetConditionGeometrically(const int i,const double max,const double min){
    const double interval = powf(max/min, 1.0/(double)(size[i]-1));
    for(int i=0;i<nreplica;i++){
    }
  }
  //*/

  template <class TResult,
	    class TFuncExchangeProb>
  void ExchangeWithNeighbor(const TResult *result,
			   TFuncExchangeProb ExchangeProb){
    assert(isInitialized);

    for(int i=0;i<nreplica;i++){
      isExchanged[i] = false;
      if(!isEven && (i/interval[dim])%size[dim]==0){
	isExchanged[i] = true;
      }
    }

    for(int i=0;i<nreplica;i++){
      int j = i + interval[dim];
      if(j>=nreplica) j -= nreplica;

      int& im = index[i]; //i(m)
      int& jm = index[j]; //j(m)
      if(isExchanged[i] || isExchanged[j]) continue;
      if(ExchangeProb(result[im],cond[i],result[jm],cond[j]) > random_number()){
	swap(cond[i],cond[j]);
	swap(im,jm);
	accept_ratio[i+dim*nreplica].Accept();
      }else{
	accept_ratio[i+dim*nreplica].Refuse();
      }
      isExchanged[i] = isExchanged[j] = true;
    }
    dim++; if(dim==NDIM) dim = 0;
    if(dim==0) isEven = !isEven;
  }

  const TCondition& GetCondition(const int i) const {return cond[i];}

  template <class TResult>
  void Output(const TResult *result){
    assert(isInitialized);
    if(output == nullptr){
      output = new std::fstream[nreplica];
      for(int i=0;i<nreplica;i++){
	const std::string filename = cond[i].FileName();
	output[i].open(filename, std::ios::out | std::ios::binary);
	if(output[i].fail())
	  std::cerr << "error: failed to open " << filename << std::endl;
      }
      for(int i=0;i<nreplica;i++)
	output[index[i]].write((char*)&cond[i],sizeof(TCondition));
    }

    for(int i=0;i<nreplica;i++)
      output[index[i]].write((char*)&result[i],sizeof(TResult));
  }

  void DebugPrint(){
    for(int i=0;i<nreplica;i++){
      for(int d=1;d<NDIM;d++)
	if(i%interval[d] == 0) std::cout << std::endl;
      std::cout << " " << index[i] << '(';
      for(int d=0;d<NDIM;d++){
	if(d!=0) std::cout << ',';
	std::cout << cond[i][d];
      }
      std::cout << ')';
    }
  }
};
