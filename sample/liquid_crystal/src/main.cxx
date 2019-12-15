#include "condition.hxx"

#include <vector>
#include <string>

#include <cassert>

#include <wham.hxx>

#define NDIM 2

int main(int argc, char** argv){
  WHAM<2,Condition,Result>::WHAMParam p;
  p.ncond = argc - 1;

  // interval of E and V for making histogram
  p.nbins[0] = 100;
  p.nbins[1] = 100;

  // number of OpenMP threads (currently, should be 1)
  p.nthreads = 1;

  // threshold to stop iteration of DoS-G calculation
  p.threshold = 0.000001;
  // number of plot for output physical quantity file
  p.npoint[0] = 1;
  p.npoint[1] = 100;

  p.print();

  std::vector<std::string> filenames;
  for(int i=1;i<argc;i++)
    filenames.push_back(argv[i]);

  WHAM<2,Condition,Result> wham(p);
  wham.ReadAscii(filenames);

  while(!wham.isIterationEnough()){
    wham.CalcDoS();
    wham.CalcFreeEnergy();
  }
  wham.CalcPhysicalValue("phys_value.dat");
  wham.OutputFreeEnergy();
  //wham.OutputDoS("dos_test.dat");
}
