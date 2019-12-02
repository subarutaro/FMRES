#ifndef HXX_WHAM
#define HXX_WHAM

#include <vector.hxx>
#include <array.hxx>
#include <sum_log.hxx>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>

#include <cmath>

//#include <omp.h>

template <int N,class TCondition,class TResult>
class WHAM{
public:
  class WHAMParam{
  public:
    int ncond;
    // for WHAM
    int nthreads;

    Vector<int,N> nbins;
    double threshold;
    Vector<int,N> npoint;

    Vector<double,N> maxs = 0.0;
    Vector<double,N> mins = 0.0;
    Vector<double,N> tics = 0.0;

    WHAMParam(){};
    WHAMParam(const WHAMParam& _p):
      ncond(_p.ncond),
      nbins(_p.nbins),
      nthreads(_p.nthreads),
      threshold(_p.threshold),
      npoint(_p.npoint){}

    void print(std::ostream& os = std::cout){
      os << "ncond:\t"   << ncond << std::endl;
      os << "nbins:\t"   << nbins << std::endl;
    };
  };
  WHAMParam p;

  TCondition    *cond  = nullptr;
  Array<int,N>  *hist  = nullptr;
  unsigned long *nhist = nullptr;

  TResult *time_average     = nullptr;
  double  *time_count       = nullptr;
  double  *free_energy      = nullptr;
  double  *prev_free_energy = nullptr;

  Array<double,N> dos;
  double dos_min;
  Array<TResult,N> average;
  Array<double ,N> count;
  TResult average_min;
public:
  WHAM(){}
  WHAM(const WHAMParam& _p) : p(_p){
    cond  = new TCondition[p.ncond];
    hist  = new Array<int,N>[p.ncond];
    nhist = new unsigned long[p.ncond];

    time_average     = new TResult[p.ncond];
    time_count       = new double[p.ncond];
    free_energy      = new double[p.ncond];
    prev_free_energy = new double[p.ncond];

    for(int i=0;i<p.ncond;i++){
      hist[i].Allocate(p.nbins);
      nhist[i] = 0;

      time_average[i] = 0.0;
      free_energy[i] = 0.0;
      prev_free_energy[i] = std::numeric_limits<double>::max();
    }
    dos.Allocate(p.nbins);
    average.Allocate(p.nbins);
    count.Allocate(p.nbins);
    for(int i=0;i<prod(p.nbins);i++){
      average[i] = 0.0;
      count[i] = 0.0;
    }
  }

  void Increment(const int c,const TResult &r){
    Vector<double,N> val = r.getBin();
    assert(p.mins <= val && val <= p.maxs);
    Vector<int,N> index = (val - p.mins) / p.tics;
    //std::cerr << val << " " << index << std::endl;

    hist[c](index)++;
    nhist[c]++;

    average(index) += r;
    count(index)   += 1.0;

    time_average[c] += r;
    time_count[c]   += 1.0;
  }

  void ReadAscii(const std::vector<std::string> filenames){
    //#pragma omp parallel for num_threads(p.nthreads)
    p.maxs = std::numeric_limits<double>::lowest();
    p.mins = std::numeric_limits<double>::max();
    for(int c=0;c<p.ncond;c++){
      std::ifstream ifs(filenames[c],std::ios::in);

      std::string line;
      std::getline(ifs,line);//skip conditions
      //cond[c] = ctmp;
      for(;std::getline(ifs,line);){
	std::stringstream str(line);
	TResult rtmp;
	str >> rtmp;
	p.maxs = max(p.maxs,rtmp.getBin());
	p.mins = min(p.mins,rtmp.getBin());
      }
    }
    p.tics = (p.maxs - p.mins)/p.nbins;
    std::cerr << "maxs: " << p.maxs << std::endl;
    std::cerr << "mins: " << p.mins << std::endl;
    std::cerr << "tics: " << p.tics << std::endl;

    for(int c=0;c<p.ncond;c++){
      std::ifstream ifs(filenames[c],std::ios::in);

      std::string line; std::getline(ifs,line);
      std::stringstream strs_c(line);
      TCondition ctmp; strs_c >> ctmp;
      cond[c] = ctmp;
      for(;std::getline(ifs,line);){
	std::stringstream str_r(line);
	TResult rtmp; str_r >> rtmp;
	//std::cerr << line << " " << rtmp << std::endl;
	Increment(c,rtmp);
      }
    }
    for(int c=0;c<p.ncond;c++){
      std::string filename = "hist" + std::to_string(c) + ".dat";
      std::ofstream ofs(filename);
      for(int i=0;i<prod(p.nbins);i++){
	ofs << Bin(i) << " " << hist[c][i] << std::endl;
      }
      ofs << std::endl;
    }
  }
#if 1
  void ReadBinary(const std::vector<std::string> filenames){
    std::cerr << "# ReadBinary" << std::endl;
    std::cerr << "# reading all files to set parameters" << std::endl;
    p.maxs = std::numeric_limits<double>::lowest();
    p.mins = std::numeric_limits<double>::max();
    for(int c=0;c<p.ncond;c++){
      std::ifstream ifs(filenames[c],std::ios::in | std::ios::binary);
      std::cerr << "# --- reading " << filenames[c] << " ---" << std::endl;

      TCondition ctmp;
      ifs.read((char*)&cond[c],sizeof(TCondition));
#ifdef OLD_FORMAT
      int nelem;
      ifs.read((char*)&nelem,sizeof(int));
#endif
      const size_t begin = ifs.tellg();
      ifs.seekg(0,ifs.end);
      const size_t end   = ifs.tellg();
      ifs.seekg(0,ifs.beg);
      ifs.read((char*)&ctmp,sizeof(TCondition));
#ifdef OLD_FORMAT
      ifs.read((char*)&nelem,sizeof(int));
      const size_t length = (end - begin) / (sizeof(double)*NDATA_RESULT);
      TResult *buffer = new TResult[length];
      for(int l=0;l<length;l++)
	for(int i=0;i<NDATA_RESULT;i++)
	  ifs.read((char*)&buffer[l].data[i],sizeof(double));
#else
      const size_t length = (end - begin + sizeof(TResult) - 1) / sizeof(TResult);
      TResult *buffer = new TResult[length];
      ifs.read((char*)buffer,length*sizeof(TResult));
#endif
#if 0
      for(int i=0;i<3;i++){
	std::cerr << ctmp;
	std::cerr << buffer[i] << std::endl;
      }
#endif
      for(int i=0;i<length;i++){
	p.maxs = max(p.maxs,buffer[i].getBin());
	p.mins = min(p.mins,buffer[i].getBin());
      }
      delete[] buffer;
    }
    p.tics = (p.maxs - p.mins) / p.nbins;

    std::cerr << "# reading all files to struct histograms" << std::endl;
    //#pragma omp parallel for num_threads(p.nthreads)
    for(int c=0;c<p.ncond;c++){
      std::ifstream ifs(filenames[c],std::ios::in | std::ios::binary);
      std::cerr << "# --- reading " << filenames[c] << " ---" << std::endl;

      TCondition ctmp;
      ifs.read((char*)&cond[c],sizeof(TCondition));
#ifdef OLD_FORMAT
      int nelem;
      ifs.read((char*)&nelem,sizeof(int));
#endif
      const size_t begin = ifs.tellg();
      ifs.seekg(0,ifs.end);
      const size_t end   = ifs.tellg();
      ifs.seekg(0,ifs.beg);
      ifs.read((char*)&ctmp,sizeof(TCondition));
#ifdef OLD_FORMAT
      ifs.read((char*)&nelem,sizeof(int));
      const size_t length = (end - begin) / (sizeof(double)*NDATA_RESULT);
      TResult *buffer = new TResult[length];
      for(int l=0;l<length;l++)
	for(int i=0;i<NDATA_RESULT;i++)
	  ifs.read((char*)&buffer[l].data[i],sizeof(double));
#else
      const size_t length = (end - begin + sizeof(TResult) - 1) / sizeof(TResult);
      TResult *buffer = new TResult[length];
      ifs.read((char*)buffer,length*sizeof(TResult));
#endif

      for(int i=0;i<length;i++) Increment(c,buffer[i]);
      delete[] buffer;
    }
  }
#endif
  Vector<double,N> Bin(const int index){
    const Vector<double,N> offsets = cumprod(p.nbins) / p.nbins;
    Vector<double,N> bin;
    int prev = 0;
    for(int i=N-1;i>=0;i--){
      bin[i] = p.mins[i] + ((index - prev) / offsets[i]) * p.tics[i];
      prev = ((index - prev) / offsets[i]);
      prev *= offsets[i];
    }
    return bin;
  }

  void CalcDoS(){
    //std::cerr << "# --- CalcDoS ---" << std::endl;
    double dos_min_tmp[p.nthreads];
    for(int i=0;i<p.nthreads;i++) dos_min_tmp[i] = std::numeric_limits<double>::max();

    //#pragma omp parallel for num_threads(p.nthreads)
    for(int i=0;i<prod(p.nbins);i++){
      double nume = 0.0;
      double deno = 0.0;
      for(int c=0;c<p.ncond;c++){
	if(hist[c][i] > 0){
	  TCondition ctmp = cond[c];
	  const Vector<double,N> bin = Bin(i);
	  const double bf = ctmp.calcLogBoltzmannFactor(bin); // returns (-\beta*H)
	  nume = sum_log(nume,log((double)hist[c][i]));
	  deno = sum_log(deno,log((double)nhist[c]) + bf + free_energy[c]);
	}
      }
      if(nume == 0.0) dos[i] = 0.0;
      else            dos[i] = nume - deno;
      if(dos[i] != 0.0) dos_min_tmp[i%p.nthreads] = (dos[i] < dos_min_tmp[i%p.nthreads]) ? dos[i]:dos_min_tmp[i%p.nthreads];
    }
    dos_min = dos_min_tmp[0];
    for(int i=1;i<p.nthreads;i++)
      dos_min = (dos_min_tmp[i] < dos_min) ? dos_min_tmp[i] : dos_min;
    /*
    for(int i=0;i<prod(p.nbins);i++)
      if(dos[i]!=0.0) std::cout << Bin(i) << " " << dos[i] << std::endl;
    //*/
    for(int i=0;i<prod(p.nbins);i++){
      if(dos[i] != 0.0) dos[i] -= dos_min;
    }
  }

  void CalcFreeEnergy(){
    //std::cerr << "# --- CalcFreeEnergy ---" << std::endl;
    for(int c=0;c<p.ncond;c++)
      prev_free_energy[c] = free_energy[c];

#pragma omp parallel for num_threads(p.nthreads)
    for(int c=0;c<p.ncond;c++){
      const Condition ctmp = cond[c];
      double tmp = 0.0;
      for(int i=0;i<prod(p.nbins);i++){
	if(dos[i]!=0.0){
	  const Vector<double,N> bin = Bin(i);
	  const double bf = ctmp.calcLogBoltzmannFactor(bin); // return (-\beta * H)
	  tmp = sum_log(tmp,dos[i] + bf);
	}
      }
      free_energy[c] = -tmp;
    }
    const double fe_min = free_energy[0];
    for(int c=0;c<p.ncond;c++) free_energy[c] -= fe_min;

    //for(int c=0;c<p.ncond;c++) std::cout << ' ' << fabs(free_energy[c] - prev_free_energy[c]); std::cout << std::endl;
    //for(int c=0;c<p.ncond;c++) std::cout << std::setprecision(12) << free_energy[c] << std::endl;

  }

  void CalcPhysicalValue(const std::string filename){
    //std::cerr << "# --- CalcPhysicalValue ---" << std::endl;
    std::ofstream ofs(filename,std::ios::out);

    average_min = std::numeric_limits<float>::max();
    for(int i=0;i<prod(p.nbins);i++){
      if(count[i] > 0.0){
	average[i] = average[i] / count[i];
	average_min = min(average_min,average[i]);
	//std::cerr << average[i] << " " << average_min << std::endl;
      }
    }
    Condition cmax,cmin;
    for(int i=0;i<N;i++){
      cmax[i] = std::numeric_limits<double>::lowest();
      cmin[i] = std::numeric_limits<double>::max();
    }
    for(int c=0;c<p.ncond;c++){
      for(int i=0;i<N;i++){
	cmax[i] = std::max(cmax[i],cond[c][i]);
	cmin[i] = std::min(cmin[i],cond[c][i]);
      }
    }
    Vector<double,N> ctics;
    for(int i=0;i<N;i++) ctics[i] =  (cmax[i] - cmin[i]) / p.npoint[i];
    Vector<int,N> coffsets = cumprod(p.npoint)/p.npoint;
    for(int c=0;c<prod(p.npoint);c++){
      int prev = 0;
      TCondition ctmp;
      for(int i=N-1;i>=0;i--){
	ctmp[i] = cmin[i] + ((c - prev) / coffsets[i]) * ctics[i];
	prev = ((c - prev) / coffsets[i]) * coffsets[i];
      }
      double pf = 0.0;
      TResult ave,sq,p4;
      ave = sq = p4 = 0.0;
      //#pragma omp declare reduction(sum_log : TResult : omp_out=sum_log(omp_out,omp_in))
      //#pragma omp declare reduction(sum_log : double  : omp_out=sum_log(omp_out,omp_in))
      //#pragma omp parallel for reduction(sum_log:ave,sq,p4,pf) num_threads(p.nthreads)
      for(int i=0;i<prod(p.nbins);i++){
	if(dos[i]!=0.0){
	  assert(count[i] > 0.0);
	  const Vector<double,N> bin = Bin(i);
	  const double bf = ctmp.calcLogBoltzmannFactor(bin);
	  const TResult a = average[i];
	  const TResult a_abs = average[i] - average_min + std::numeric_limits<float>::min();

	  ave = sum_log(ave,log(a_abs)   + dos[i] + bf);
	  sq  = sum_log(sq, log(a*a)     + dos[i] + bf);
	  p4  = sum_log(p4, log(a*a*a*a) + dos[i] + bf);
	  pf  = sum_log(pf,                dos[i] + bf);
	}
	//std::cout << ave << " " << pf << std::endl;
      }
      ave = exp(ave - pf) + average_min;
      sq  = exp(sq - pf);
      p4  = exp(p4 - pf);
      //std::cout << " " << ctmp << " " << ave << " " << average_min << " " << pf << std::endl;

      for(int d=1;d<N;d++) if(c%coffsets[d]==0) ofs << std::endl;
      ofs << " " << ctmp
	  << " " << ave
	//<< " " << sq
	//<< " " << p4
	  << " " << pf
	  << std::endl;
    }
  }

  bool isIterationEnough(){
    double error = std::numeric_limits<double>::lowest();
    double fe_max = 0.0;
    for(int i=1;i<p.ncond;i++){ // free_energy[0] must be 0
      //const double diff = fabs((free_energy[i] - prev_free_energy[i])/free_energy[i]);
      const double diff = fabs(free_energy[i] - prev_free_energy[i]);
      error = std::max(error,diff);
      fe_max = std::max(fe_max, fabs(free_energy[i]));
    }
    const double threshold = fe_max * p.threshold;
    //std::cerr << "# max error is " << error << " > " << threshold << std::endl;
    return (error < threshold);
  }

  void OutputHistogram(const std::string prefix){
    for(int c=0;c<p.ncond;c++){
      std::stringstream strs;
      strs << "hist" << std::setw(6) << std::setfill('0') << c << ".dat";
      const std::string filename = strs.str();
      std::ofstream ofs(filename);
      assert(!ofs.fail());
      ofs << "#" << cond[c] << std::endl;
      for(int i=0;i<prod(p.nbins);i++){
	ofs << Bin(i) << " " << hist[c][i] << std::endl;
      }
    }
  }

  void OutputFreeEnergy(const std::string prefix = "./"){
    std::string filename = prefix + "free_energy.dat";
    std::ofstream ofs(filename);
    for(int c=0;c<p.ncond;c++){
      ofs << cond[c] << " " << std::setprecision(15) << free_energy[c] << std::endl;
    }
  }

  void OutputDoS(const std::string prefix = "./"){
    std::string filename = prefix + "dens_state.dat";
    std::ofstream ofs(filename);
    for(int i=0;i<prod(p.nbins);i++){
      Vector<double,N> bin = Bin(i);
      ofs << bin[i] << " " << dos[i] << std::endl;
    }
  }

  void ReadFreeEnergy(const std::string filename){
    std::ifstream ifs(filename);
    std::string line;
    int count = 0;
    for(std::string line;std::getline(ifs,line);){
      count++;
      std::cout << line << std::endl;
    }
    assert(p.ncond == count);
    ifs.clear();
    ifs.seekg(0,ifs.beg);
    for(int c=0;c<p.ncond;c++){
      std::string line;
      std::getline(ifs,line);
      std::stringstream strs(line);
      Condition dummy;
      strs >> dummy >> free_energy[c];
    }
  }
};

#endif
