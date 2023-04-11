// The CCfits headers are expected to be installed in a subdirectory of
// the include path.
// The <CCfits> header file contains all that is necessary to use both the CCfits
// library and the cfitsio library (for example, it includes fitsio.h) thus making
// all of cfitsio’s macro definitions available.
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
// this includes 12 of the CCfits headers and will support all CCfits operations.
// the installed location of the library headers is $(ROOT)/include/CCfits
// to use the library either add -I$(ROOT)/include/CCfits or #include <CCfits/CCfits>
// in the compilation target.
#include <stddef.h>
#include <limits.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <algorithm>
#include <numeric>

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifndef INFINITY
#define INFINITY (1.0 / 0.0)
#endif

// The library is enclosed in a namespace.
using namespace std;
//using namespace stapl;

//typedef std::numeric_limits< double > dbl;
const float EPSILON = 0.00000001f;
using namespace std;


void myBls( vector<double> scannedWeights ,vector<double> scannedWeightedFlux,
           vector<double> time, double helper_d, long long size){
  double r,s, d;
  int min_i1=-1, min_i2=-1;
  double d_min= DBL_MAX;
  double ex_time = DBL_MIN;

  //double mem_time = 0.0;
  //double global_mem_time = 0.0;
  std::setprecision(15);
  int p = atoi(getenv("OMP_NUM_THREADS"));
  cout << "found " << p << " threads in the program." << endl;
  long long  totalArea = (size*size)/2;
  long long  areaPerProc = totalArea/p;
  cout << "computed " << areaPerProc << " areaPerProc." << endl;
  vector<long long> portion;
  portion.push_back(0);
  portion.push_back(sqrt(2*areaPerProc));
  long long areaCovered = (portion[1]*portion[1])/2;
  for(size_t i=1; i< (size_t) p-1; i++){
      portion.push_back(sqrt(2*(areaPerProc+areaCovered)));
      areaCovered = (portion[i+1]*portion[i+1])/2;
  }
  portion.push_back(size);
   std::reverse(portion.begin(),portion.end());
  for(size_t i=0; i< (size_t) p+1; i++){
      portion[i] = size - portion[i] ;
      cout << "computed area offset for " << i << ":\t" << portion[i] << endl;
  }
  #pragma omp parallel  shared(scannedWeights, scannedWeightedFlux) private(r,s,d) reduction(min:d_min) reduction(max:ex_time )
  {
      double wtime = omp_get_wtime();
    int num_loc = omp_get_num_threads();
    int loc_id = omp_get_thread_num();
    double reg1,reg2;
    //long long lower_limit = (loc_id==0) ? 0 : portion[loc_id-1];
    for(long long i1=portion[loc_id]; i1< portion[loc_id+1]; i1++){
        //double st_mem_time = omp_get_wtime();
        reg1= scannedWeights[i1];
        reg2= scannedWeightedFlux[i1];
        //double after_mem_time = omp_get_wtime();
        //mem_time+= after_mem_time - st_mem_time;
        for(long long i2=(i1+1); i2< (long long) size; i2++){
            //st_mem_time = omp_get_wtime();
            r = scannedWeights[i2] - reg1;
            s = scannedWeightedFlux[i2] - reg2;
            //r = reg1+ 2.0*reg2 - reg1;
            //s = 4.0*reg2 - reg2;
            //after_mem_time = omp_get_wtime();
            //mem_time+= after_mem_time - st_mem_time;

            d = (helper_d - ( (s*s) / ( 1.0 * r *  (1.0-r) )));
            //double redctime = omp_get_wtime();
            if(d < d_min){
                d_min = d;
                min_i1=i1;
                min_i2=i2;
            }
        }
    }
    wtime = omp_get_wtime() - wtime;
    // reduce global times
    if (wtime > ex_time ){
        ex_time=wtime; //global_mem_time=mem_time;
    }
  }
  if(min_i1!=-1 && min_i2!=-1) {
      double corrected_p =  -0.002*(time[min_i2] - time[min_i1]) + 0.032;
    cout << std::fixed << "resulting i1: " << min_i1 << "\ti2: " << min_i2 << "\td: " << d_min  << " time: " << ex_time << " period: "  <<  corrected_p << '\n' ;
    //cout << std::fixed << " mem-time: " << global_mem_time   << '\n' ;
  }
  else cout << "could not find any pairs. latest d: " << d << endl;

}
int main(int argc, char **argv)
{
 // cout.precision(dbl::max_digits10);

  long long size = atoi(argv[1]);
  vector<double> view_flux;
  vector<double> view_fluxerr;
  vector<double> view_time;
  for(size_t i=0; i< (size_t) size; i++){
    double flux =  (double)rand() / RAND_MAX;
    view_flux.push_back(0.005 + flux * (3.2 - 0.005));
    view_fluxerr.push_back(0.0003);
    view_time.push_back(view_flux[i] + 0.002);
  }
  cout << "created inputs\t" << endl << flush;


  vector<double> view_squaredfluxerr;
  // square errors
  transform(view_fluxerr.begin(), view_fluxerr.end(),  std::back_inserter(view_squaredfluxerr),
            [](double x) {  return pow((double)x,-2.0); });
  // sum squared errors and take power of sum
  double sumW = pow(std::accumulate(view_squaredfluxerr.begin(), view_squaredfluxerr.end(),  0), -1);
  cout << "computed sumW \t" << sumW << endl << flush;

  vector<double> view_weight;
  // create weight list with sumw*err^2
  transform(view_squaredfluxerr.begin(), view_squaredfluxerr.end(),  std::back_inserter(view_weight),
            [sumW](double x) {   return sumW*x; });

  vector<double> view_weightedFlux;
  transform(view_flux.begin(), view_flux.end(), view_weight.begin(),  std::back_inserter(view_weightedFlux),
       [](double f, double w) {  return f*w; });

  vector<double>  view_d;
  transform(view_squaredfluxerr.begin(), view_squaredfluxerr.end(), view_weight.begin(), std::back_inserter(view_d),
            [](double sf, double w) {   return w * pow(sf,2.0); });
  double helper_d = accumulate(view_d.begin(), view_d.end(), 0);

  // prefix sum weighted flux
  vector<double> v_scanned_weightedFlux;
  v_scanned_weightedFlux.push_back(view_weightedFlux[0]);
  for(size_t i=1; i< (size_t) view_weightedFlux.size(); i++){
    v_scanned_weightedFlux.push_back(view_weightedFlux[i] + v_scanned_weightedFlux[i-1]);
  }
  // prefix sum weights
  vector<double> v_scanned_weights;
  v_scanned_weights.push_back(view_weight[0]);
  for(size_t i=1; i< (size_t) view_weight.size(); i++){
    v_scanned_weights.push_back(view_weight[i] + v_scanned_weights[i-1]);
  }

  //double start_time = omp_get_wtime();
  // do once. create time, flux and flux error.
  myBls(v_scanned_weights, v_scanned_weightedFlux, view_time,  helper_d,  size);
  //double end_time = omp_get_wtime();

  //std::cout <<  "-finished BLS with time : " << end_time-start_time << "\n";
  //double max_time=0.0;
  //MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  //stapl::do_once(
  //        [&]{ std::cout << "\nmax time : " << max_time << std::endl;}
  //);

  return 0;
}
