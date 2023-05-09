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
#include <CCfits/CCfits>
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
using namespace CCfits;
using namespace std;
//using namespace stapl;

typedef std::numeric_limits< double > dbl;
const float EPSILON = 0.00000001f;
using namespace std;

 size_t preprocess(std::string fileName,  vector<double> &fluxBLS, vector<double> &fluxErrBLS,  vector<double> &timeBLS)
{
    // read a table and explicitly read selected columns. To read instead all the
    // data on construction, set the last argument of the FITS constructor
    // call to ’true’. This functionality was tested in the last release.

    std::vector<string> hdus(2);
    hdus[0] = "PRIMARY";
    hdus[1] = "LIGHTCURVE";
    //hdus[2] = "APERTURE";
    std::auto_ptr<FITS> pInfile(new FITS(fileName,Read,hdus,true));
    ExtHDU& table = pInfile->extension(hdus[1]);
    size_t fileSize = table.axis(1);
    cout << "read filesize: " << fileSize << endl;
    std::valarray <double> time, flux, flux_err;
    table.column("TIME").read( time, 1,fileSize );
    table.column("FLUX").read( flux, 1,fileSize );
    table.column("FLUX_ERR").read( flux_err, 1,fileSize );

    cout.precision(dbl::max_digits10);

    size_t toBeDeleted=0;
    for(size_t i=0; i<time.size(); i++){

        if(std::isnan(time[i]) || std::isnan(flux[i]) || std::isnan(flux_err[i]) ){
            toBeDeleted++;
        }
    }
    cout << "cleaning NaNs: " << toBeDeleted << endl;

    size_t newfileSize = time.size() - toBeDeleted;
    vector<double> tempflux, temperr, temptime;
    for(size_t i=0; i< (size_t) time.size(); i++){
        if(std::isnan(time[i]) || std::isnan(flux[i]) || std::isnan(flux_err[i])) continue;
        tempflux.push_back(flux[i]);
        temperr.push_back(flux_err[i]);
        temptime.push_back(time[i]);
    }

    /*
    cout << "before normalizing: TIME \t FLUX \t FLUXERR" << endl;
    for(size_t i=0; i< (size_t) newfileSize; i++){
        cout << temptime[i] << '\t' << tempflux[i] << '\t' <<  temperr[i] << '\n';
    }
     */

    // for flux
    double min_err, max_err, min_fl, max_fl;
    min_fl =  *std::min_element(tempflux.begin(),tempflux.end());
    max_fl =  *std::max_element(tempflux.begin(),tempflux.end());
    min_err =  *std::min_element(temperr.begin(),temperr.end());
    max_err =  *std::max_element(temperr.begin(),temperr.end());
    std::vector<double> normFlux(newfileSize), normFluxErr(newfileSize);
    std::transform(tempflux.begin(), tempflux.end(), normFlux.begin(),
                   [min_fl, max_fl](double x) { return (x - min_fl) / (max_fl-min_fl); });

    std::transform(temperr.begin(), temperr.end(), normFluxErr.begin(),
                   [min_err, max_err](double x) { return (x - min_err) / (max_err-min_err); });

    // print after norm
    /*
    cout << "AFTER NORM: FLUX \t FLUXERR" << endl;
    for(size_t i=0; i< (size_t) newfileSize; i++){
        cout << std::fixed << normFlux[i] << '\t' <<  normFluxErr[i] << '\n';
    }
     */

    size_t numzeros=0;
    // clean 0.0s so that power does not crash the program with infs
     for(size_t i=0; i< (size_t) newfileSize; i++){
        if(((normFlux[i]<EPSILON) && (normFlux[i]>= -EPSILON)) || ((normFluxErr[i]<EPSILON) && (normFluxErr[i]>= -EPSILON)) ){
            numzeros++;
        }
    }
     cout << "cleaning 0.0s: " << numzeros << endl;
     newfileSize -= numzeros;

     //fluxBLS = new double [newfileSize];
     //fluxErrBLS = new double [newfileSize];
     //timeBLS = new double [newfileSize];
     //int ctr=0;
     // original size of norms contain 0s too, iterate over the previous size
      for(size_t i=0; i< (size_t) newfileSize+numzeros; i++){
            if(((normFlux[i]<EPSILON) && (normFlux[i]>= -EPSILON)) || ((normFluxErr[i]<EPSILON) && (normFluxErr[i]>= -EPSILON)) ) continue;
                fluxBLS.push_back(normFlux[i]);
                fluxErrBLS.push_back(normFluxErr[i]);
                timeBLS.push_back(temptime[i]);
            }

    /*
    cout << "TIME \t FLUX \t FLUXERR" << endl;
    for(size_t i=0; i< (size_t) newfileSize; i++){
    cout << std::fixed<< timeBLS[i] << '\t' << fluxBLS[i] << '\t' <<  fluxErrBLS[i] << '\n';
    }
    */

    return newfileSize;
}


void myBls( vector<double> scannedWeights ,vector<double> scannedWeightedFlux,
           vector<double> time, double helper_d, size_t size){
  double r,s, d;
  int min_i1=-1, min_i2=-1;
  double d_min= DBL_MAX;
  double reg1, reg2;
  //cout.precision(dbl::max_digits10);

    double wtime = omp_get_wtime();
    for(long long i1=0; i1< size; i1++){
        reg1= scannedWeights[i1];
        reg2= scannedWeightedFlux[i1];
      for(long long i2=(i1+1); i2< size; i2++){
        //double loop_time = omp_get_wtime();
        r = scannedWeights[i2] - reg1;
        s = scannedWeightedFlux[i2] - reg2;

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

    if(min_i1!=-1 && min_i2!=-1) {
        double corrected_p =  -0.002*(time[min_i2] - time[min_i1]) + 0.032;
        double p =  time[min_i2] - time[min_i1];
        cout << std::fixed << "resulting i1: " << min_i1 << "\ti2: " << min_i2 << "\td: " << d_min  << " time: " << wtime << " corrected period: "  <<  corrected_p  << " period: "  <<  p << " time[i2]: " << time[min_i2] << " time[i1]: " << time[min_i1] << '\n' ;
    }
    else cout << "could not find any pairs. latest d: " << d << endl;

}
int main(int argc, char **argv)
{
   FITS::setVerboseMode(true);
    cout.precision(dbl::max_digits10);
    std::string fileName = argv[1];

    size_t size ;
    vector<double> view_flux;
    vector<double> view_fluxerr;
    vector<double> view_time;

     size = preprocess(fileName, view_flux, view_fluxerr, view_time);
    //std::cout  << ": successfully readImage() with size: " << size << std::endl;

    vector<double> view_squaredfluxerr;
    // square errors
    transform(view_fluxerr.begin(), view_fluxerr.end(),  std::back_inserter(view_squaredfluxerr),
                   [](double x) { cout.precision(dbl::max_digits10); return 1/(x *x ); });
    // sum squared errors and take power of sum
    double sumW = 1.0/(accumulate(view_squaredfluxerr.begin(), view_squaredfluxerr.end(),  0));
    cout << "computed sumW \t" << sumW << endl << flush;

    vector<double> view_weight;
    // create weight list with sumw*err^2
    transform(view_squaredfluxerr.begin(), view_squaredfluxerr.end(),  std::back_inserter(view_weight),
                   [sumW](double x) {  cout.precision(dbl::max_digits10); return sumW*x; });

    // for s
    vector<double> view_weightedFlux;
    transform(view_flux.begin(), view_flux.end(), view_weight.begin(),  std::back_inserter(view_weightedFlux),
                   [](double f, double w) {  cout.precision(dbl::max_digits10); return f*w; });

    // for d
    vector<double>  view_d;
    transform(view_flux.begin(), view_flux.end(), view_weight.begin(), std::back_inserter(view_d),
                   [](double f, double w) {  cout.precision(dbl::max_digits10); return w * f * f; });
    double helper_d = accumulate(view_d.begin(), view_d.end(), 0);
    cout << "constant sum for d : " << helper_d << endl;

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
