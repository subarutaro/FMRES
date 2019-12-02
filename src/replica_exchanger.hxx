#include <cmath>
#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <cstring>

#include <vector.hxx>
#include <array.hxx>

//#include "permutation.hxx"

double random_number(){
  return (double)rand() / ((double)RAND_MAX + 1.0);
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
  Vector<int,NDIM> size;
  Vector<int,NDIM> offset;

  int nreplica;
  bool isInitialized = false;
  Array<AcceptRatio,NDIM> accept_ratio[NDIM];

  int dim = 0;
  bool isEven = true;

  std::fstream *output = nullptr;
  std::fstream of_mtoi;
  std::fstream of_itom;
public:
  Array<int,NDIM> mtoi;
  Array<int,NDIM> itom;
  Array<TCondition,NDIM> cond;
  Array<TCondition,NDIM> next;

  ReplicaExchanger(){
    size = 1;
    nreplica = 1;
  };
  ~ReplicaExchanger(){
    if(of_itom.is_open()) of_itom.close();
    if(of_mtoi.is_open()) of_mtoi.close();
    /*
    if(output != nullptr)
      for(int i=0;i<nreplica;i++)
	if(outout[i].is_open()) output[i].close();
    //*/
  }
  void Initialize(const Vector<int,NDIM> _size){
    assert(!isInitialized);
    assert(_size>0);

    size = _size;
    nreplica = prod(size);
    offset = cumprod(size)/size;

    cond.Allocate(size);
    next.Allocate(size);
    mtoi.Allocate(size);
    itom.Allocate(size);
    for(int i=0;i<nreplica;i++) mtoi[i] = i;
    for(int i=0;i<nreplica;i++) itom[i] = i;
    for(int d=0;d<NDIM;d++) accept_ratio[d].Allocate(size);

    dim = 0;
    isEven = true;
    isInitialized = true;
  };

  void SetConditionRegularInterval(const Vector<double,NDIM> max,const Vector<double,NDIM> min){
    Vector<double,NDIM> diff;
    for(int i=0;i<NDIM;i++) diff[i] = (size[i]>1)?(max[i] - min[i]) / (double)(size[i]-1) : 0.0;

    for(int i=0;i<nreplica;i++){
      Vector<int,NDIM> indices = Index2Indices(i);
      for(int d=0;d<NDIM;d++) cond[i][d] = min[d] + diff[d]*indices[d];
    }
  }
  //*
    void SetConditionGeometrically(const int i,const double max,const double min){
    const double interval = powf(max/min, 1.0/(double)(size[i]-1));
    for(int i=0;i<nreplica;i++){
    }
  }
  //*/

  bool OptimizeConditionsOfDimension(const int  _dim,
				     const double ideal_ratio = 0.2){
    double max = 0.0, min = 1.0;
    for(int rep=0;rep<nreplica;rep++){
      const Vector<int,NDIM> i = Index2Indices(rep);
      if(i[_dim] == size[_dim]-1) continue;
      const double ar = accept_ratio[_dim](i).Get();
      if(ar > max) max = ar;
      if(ar < min) min = ar;
    }
    if(max < ideal_ratio){
      std::cerr << "error: max accept ratio(=" << max << ") is smaller than ideal ratio (=" << ideal_ratio << ").";
      exit(EXIT_FAILURE);
    }
    if(min > ideal_ratio) return true;

    for(int rep=0;rep<nreplica;rep++){
      const Vector<int,NDIM> i = Index2Indices(rep);
      if(i[_dim] == size[_dim]-1) continue;
      Vector<int,NDIM> j = i;
      j[_dim] += 1;
      const double ar = accept_ratio[dim](i).Get();
      const double dc = (cond(i)[_dim] - cond(j)[_dim]);
      if(ar < ideal_ratio){
	if(cond(i)[_dim] != 0 && cond(i)[_dim] != size[_dim]-1) cond(i)[_dim] += dc*0.01;
	if(cond(j)[_dim] != 0 && cond(j)[_dim] != size[_dim]-1) cond(j)[_dim] -= dc*0.01;
      }else{
	if(cond(i)[_dim] != 0 && cond(i)[_dim] != size[_dim]-1) cond(i)[_dim] -= dc*0.01;
	if(cond(j)[_dim] != 0 && cond(j)[_dim] != size[_dim]-1) cond(j)[_dim] += dc*0.01;
      }
    }
    return false;
  }

  Vector<int,NDIM> Index2Indices(const int index){
    Vector<int,NDIM> indices;
    int prev = 0;
    for(int d=NDIM-1;d>=0;d--){
      indices[d] = (index-prev)/offset[d];
      prev += indices[d] * offset[d];
    }
    return indices;
  }

  int Indices2Index(Vector<int,NDIM> indices){
    return indices % offset;
  }

  template <class TResult,
	    class TFuncExchangeProb>
    void ExchangeWithNeighbor(const TResult *result,
			      TFuncExchangeProb ExchangeProb){
    assert(isInitialized);
    if(nreplica == 1) return;

    static unsigned long count = 0;

    for(int i=0;i<nreplica;i++){
      if(size[dim]<2) break;

      Vector<int,NDIM> ivec = Index2Indices(i);
      if(((ivec[dim]%2)!=(count%2))) continue;

      Vector<int,NDIM> jvec = ivec;
      jvec[dim]++;
      //if(jvec[dim] == size[dim]) jvec[dim] = 0;
      if(jvec[dim] == size[dim]) continue;

      int& im = mtoi(ivec); //i(m)
      int& jm = mtoi(jvec); //j(m)
      if(ExchangeProb(result[im],cond[im],
		      result[jm],cond[jm]) > random_number()){
	/*
	std::cout << "dim " << dim << ":"
		  << "(" << ivec << ":" << cond[im] << ")"
		  << " <--> "
		  << "(" << jvec << ":" << cond[jm] << ")"
		  << std::endl;
	//*/
	std::swap(cond[im],cond[jm]);
	std::swap(im,jm);
	accept_ratio[dim](ivec).Accept();
      }else{
	accept_ratio[dim](ivec).Refuse();
      }
    }
    for(int m=0;m<nreplica;m++) itom[mtoi[m]] = m;

    dim++;
    if(dim==NDIM){
      dim -= NDIM;
      count++;
    }

    for(int i=0;i<nreplica;i++){
      if(size[dim]<2) break;
      const Vector<int,NDIM> ivec = Index2Indices(i);
      if((ivec[dim]%2)!=(count%2)) continue;
      Vector<int,NDIM> jvec = ivec;
      jvec[dim]++;
      if(jvec[dim]==size[dim]) jvec[dim] = 0;

      next[mtoi(ivec)] = cond[mtoi(jvec)];
      next[mtoi(jvec)] = cond[mtoi(ivec)];
    }
  }

  void OutputIndex(std::string prefix = ""){
    assert(isInitialized);
    static bool isInitialized = false;
    static unsigned long count = 0;

    const std::string itom_file = prefix + "itom.dat";
    const std::string mtoi_file = prefix + "mtoi.dat";

    if(!isInitialized){
      std::cerr << "# initializing itom and mtoi output ..." << std::endl;
      of_itom.open(itom_file,std::ios::out);
      if(of_itom.fail()){
	std::cerr << "error: failed to open " << itom_file << std::endl;
	exit(EXIT_FAILURE);
      }
      for(int d=0;d<NDIM;d++)
	of_itom << "# size of dim " << d << " = " << size[d] << std::endl;

      of_mtoi.open(mtoi_file,std::ios::out);
      if(of_mtoi.fail()){
	std::cerr << "error: failed to open " << mtoi_file << std::endl;
	exit(EXIT_FAILURE);
      }
      for(int d=0;d<NDIM;d++)
	of_mtoi << "# size of dim " << d << " = " << size[d] << std::endl;
      isInitialized = true;
    }

    of_itom << count;
    of_mtoi << count;
    for(int i=0;i<nreplica;i++){
      of_itom << " " << std::setw(5) << itom[i];
      of_mtoi << " " << std::setw(5) << mtoi[i];
    }
    of_itom << std::endl;
    of_mtoi << std::endl;
    count++;
  }

  template <class TResult>
  void OutputResult(const TResult *result,const std::string prefix = "",int nrep = 0,const int offset = 0){
    if(nrep == 0) nrep = nreplica;
    if(output == nullptr){
      output = new std::fstream[nrep];
      for(int m=0;m<nrep;m++){
	std::cout << "# Initializing " << cond[mtoi[offset+m]].FileName() << " ..." << std::endl;
	output[m].open(prefix+cond[mtoi[offset+m]].FileName(), std::ios::out | std ::ios::binary);
	output[m].write((char*)&cond[mtoi[offset+m]],sizeof(TCondition));
      }
    }
    #pragma omp parallel for
    for(int m=0;m<nrep;m++)
      output[m].write((char*)&result[mtoi[offset+m]],sizeof(TResult));
  }

  void OutputAcceptanceRatio(std::string prefix = ""){
    const std::string filename = prefix+"acceptance_ratio.dat";
    std::ofstream ofs(filename);
    if(ofs.fail()){
      std::cerr << "error: open " << filename << " failed" << std::endl;
      exit(EXIT_FAILURE);
    }

    for(int dim=0;dim<NDIM;dim++){
      for(int i=0;i<nreplica;i++){
	Vector<int,NDIM> ivec = Index2Indices(i);
	Vector<int,NDIM> jvec = ivec;
	jvec[dim]++;
	if(jvec[dim]==size[dim]) jvec[dim] = 0;
	const TCondition ci = cond[mtoi(ivec)];
	const TCondition cj = cond[mtoi(jvec)];
	for(int d=0;d<NDIM;d++) ofs << " " << ci[d];
	ofs << "  ";
	for(int d=0;d<NDIM;d++) ofs << " " << cj[d];
	ofs << "  " << accept_ratio[dim](ivec).Get();

	for(int d=0;d<NDIM;d++)
	  if(i%offset[d] == offset[d]-1) ofs << std::endl;
      }
    }
  }

  void DebugPrint(){
    for(int i=0;i<nreplica;i++){
      for(int d=1;d<NDIM;d++)
	if(i%offset[d] == 0) std::cout << std::endl;
      std::cout << " " << mtoi[i] << '(';
      for(int d=0;d<NDIM;d++){
	if(d!=0) std::cout << ',';
	std::cout << cond[i][d];
      }
      std::cout << ')';
    }
  }
};
