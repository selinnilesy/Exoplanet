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
#include <papi.h>

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

  typedef struct {
    double d;
    long long min_i1;
    long long min_i2;
  } solution_t;

solution_t solution = { DBL_MAX, -1,-1 };
#pragma omp declare reduction(get_mind : solution_t :\
    omp_out = omp_out.d < omp_in.d ? omp_out : omp_in)\
    initializer (omp_priv=(omp_orig))

  //double mem_time = 0.0;
  //double global_mem_time = 0.0;
  std::setprecision(15);

    int retval, EventSet = PAPI_NULL;
long long values [1];

/* Initialize the PAPI library */
retval = PAPI_library_init(PAPI_VER_CURRENT);
if(retval != PAPI_VER_CURRENT)
{   fprintf (stderr, "PAPI library init error!\n");
    exit();
}
/* Create the Event Set */
if (PAPI_create_eventset(&EventSet) != PAPI_OK)   handle_error(1);
/* Add Total Instructions Executed to our EventSet */
if (PAPI_add_event(EventSet, PAPI_L1_DCM) != PAPI_OK) handle_error (1);
/* Start counting */
if (PAPI_start(EventSet) != PAPI_OK) handle_error(1);

  double wtime = omp_get_wtime();
  #pragma omp parallel  shared(scannedWeights, scannedWeightedFlux) private(r,s,d) reduction(get_mind:solution )
  {
    int num_loc = omp_get_num_threads();
    int loc_id = omp_get_thread_num();
    double reg1,reg2;
    #pragma omp for schedule(dynamic,1)
    for(long long i1=0; i1< size; i1++){
        reg1= scannedWeights[i1];
        reg2= scannedWeightedFlux[i1];
        for(long long i2=(i1+1); i2< size; i2++){
            r = scannedWeights[i2] - reg1;
            s = scannedWeightedFlux[i2] - reg2;

            d = (helper_d - ( (s*s) / ( 1.0 * r *  (1.0-r) )));
            solution_t temp_solution = { d, i1,i2 };
           if(temp_solution.d < solution.d){
                    solution.d = temp_solution.d;
                    solution.min_i1=temp_solution.min_i1;
                    solution.min_i2=temp_solution.min_i2;
           }
        }
    }
  }
  ex_time = omp_get_wtime() - wtime;

    if (PAPI_read(EventSet, values) != PAPI_OK) handle_error(1);
    if (PAPI_stop(EventSet, values) != PAPI_OK) handle_error(1);

  if(solution.min_i1!=-1 && solution.min_i2!=-1) {
      double corrected_p =  -0.002*(time[solution.min_i2] - time[solution.min_i1]) + 0.032;
    cout << std::fixed << "resulting i1: " << solution.min_i1 << "\ti2: " << solution.min_i2 << "\td: " << solution.d  << " time: " << ex_time << " period: "  <<  corrected_p << '\n' ;
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
  for(long long i=0; i< size; i++){
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
            [](double sf, double w) {   return w * (1.0 / sf); });
  double helper_d = accumulate(view_d.begin(), view_d.end(), 0);

  // prefix sum weighted flux
  vector<double> v_scanned_weightedFlux;
  v_scanned_weightedFlux.push_back(view_weightedFlux[0]);
  for(long long i=1; i< (long long) view_weightedFlux.size(); i++){
    v_scanned_weightedFlux.push_back(view_weightedFlux[i] + v_scanned_weightedFlux[i-1]);
  }
  // prefix sum weights
  vector<double> v_scanned_weights;
  v_scanned_weights.push_back(view_weight[0]);
  for(long long i=1; i< (long long) view_weight.size(); i++){
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
