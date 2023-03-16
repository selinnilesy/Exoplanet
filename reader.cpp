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
#include <cmath>
#include <iostream>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stapl/runtime/runtime.hpp>
#include "benchmarks/algorithms/utilities.hpp"
#include "benchmarks/algorithms/timer.hpp"
#include <algorithm>

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
const float EPSILON = 0.00000000001f;

double *timeBLS, *fluxBLS, *fluxErrBLS;

 int readTable(std::string fileName)
{
    // read a table and explicitly read selected columns. To read instead all the
    // data on construction, set the last argument of the FITS constructor
    // call to ’true’. This functionality was tested in the last release.

    std::vector<string> hdus(3);
    hdus[0] = "PRIMARY";
    hdus[1] = "LIGHTCURVE";
    hdus[2] = "APERTURE";
    std::auto_ptr<FITS> pInfile(new FITS(fileName,Read,hdus,true));
    ExtHDU& table = pInfile->extension(hdus[1]);
    int fileSize = table.axis(1);
    cout << "read filesize: " << fileSize << endl;
    std::valarray <double> time, flux, flux_err;
    table.column("TIME").read( time, 1,fileSize );
    table.column("SAP_FLUX").read( flux, 1,fileSize );
    table.column("SAP_FLUX_ERR").read( flux_err, 1,fileSize );

    cout.precision(dbl::max_digits10);

    int toBeDeleted=0;
    for(size_t i=0; i<time.size(); i++){

        if(std::isnan(time[i]) || std::isnan(flux[i]) || std::isnan(flux_err[i]) ){
            toBeDeleted++;
        }
    }
    cout << "cleaning NaNs: " << toBeDeleted << endl;

    int newfileSize = time.size() - toBeDeleted;
    vector<double> tempflux, temperr;
    for(size_t i=0; i< (size_t) time.size(); i++){
        if(std::isnan(time[i]) || std::isnan(flux[i]) || std::isnan(flux_err[i])) continue;
        tempflux.push_back(flux[i]);
        temperr.push_back(flux_err[i]);
    }

    /*
     * cout << "before normalizing: TIME \t FLUX \t FLUXERR" << endl;
    for(size_t i=0; i< (size_t) newfileSize; i++){
        cout << timeBLS[i] << '\t' << tempflux[i] << '\t' <<  temperr[i] << '\n';
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

    int numzeros=0;
    // clean 0.0s so that power does not crash the program with infs
     for(size_t i=0; i< (size_t) newfileSize; i++){
        if(((normFlux[i]<EPSILON) && (normFlux[i]>= -EPSILON)) || ((normFluxErr[i]<EPSILON) && (normFluxErr[i]>= -EPSILON)) ){
            numzeros++;
        }
    }
     cout << "cleaning also 0.0s: " << numzeros << endl;
     newfileSize -= numzeros;
     timeBLS= new double[newfileSize];
     fluxBLS= new double[newfileSize];
     fluxErrBLS= new double[newfileSize];
     int ctr=0;
    for(size_t i=0; i< (size_t) newfileSize+numzeros; i++){
        if(((normFlux[i]<EPSILON) && (normFlux[i]>= -EPSILON)) || ((normFluxErr[i]<EPSILON) && (normFluxErr[i]>= -EPSILON)) ) continue;
        fluxBLS[ctr] = normFlux[i];
        fluxErrBLS[ctr] = normFluxErr[i];
        timeBLS[ctr] = time[i];
        ctr++;
    }

    /*
     cout << "TIME \t FLUX \t FLUXERR" << endl;
    for(size_t i=0; i< (size_t) newfileSize; i++){
        cout << std::fixed<< timeBLS[i] << '\t' << fluxBLS[i] << '\t' <<  fluxErrBLS[i] << '\n';
    }
     */

    return newfileSize;
}

/*
void compute_objective(
    double y_in,
    double y_out,
    double ivar_in,
    double ivar_out,
    int obj_flag,
    double* objective,
    double* log_likelihood,
    double* depth,
    double* depth_err,
    double* depth_snr
) {
    if (obj_flag) {
        double arg = y_out - y_in;
        *log_likelihood = 0.5*ivar_in*arg*arg;
        *objective = *log_likelihood;
    } else {
        *depth = y_out - y_in;
        *depth_err = sqrt(1.0 / ivar_in + 1.0 / ivar_out);
        *depth_snr = *depth / *depth_err;
        *objective = *depth_snr;
    }
}

int run_bls (
    // Inputs
    int N,                   // Length of the time array
    double* t,               // The list of timestamps
    double* y,               // The y measured at ``t``
    double* ivar,            // The inverse variance of the y array

    int n_periods,
    double* periods,         // The period to test in units of ``t``

    int n_durations,         // Length of the durations array
    double* durations,       // The durations to test in units of ``bin_duration``
    int oversample,          // The number of ``bin_duration`` bins in the maximum duration

    int obj_flag,            // A flag indicating the periodogram type
                             // 0 - depth signal-to-noise
                             // 1 - log likelihood

    // Outputs
    double* best_objective,  // The value of the periodogram at maximum
    double* best_depth,      // The estimated depth at maximum
    double* best_depth_err,  // The uncertainty on ``best_depth``
    double* best_duration,   // The best fitting duration in units of ``t``
    double* best_phase,      // The phase of the mid-transit time in units of
                             // ``t``
    double* best_depth_snr,  // The signal-to-noise ratio of the depth estimate
    double* best_log_like    // The log likelihood at maximum
) {
    // Start by finding the period and duration ranges
    double max_period = periods[0], min_period = periods[0];
    int k;
    for (k = 1; k < n_periods; ++k) {
        if (periods[k] < min_period) min_period = periods[k];
        if (periods[k] > max_period) max_period = periods[k];
    }
    cout << "min_period: " << min_period << endl ;
    cout << "max_period: " << max_period << endl ;
    if (min_period < DBL_EPSILON) return 1;
    double min_duration = durations[0], max_duration = durations[0];
    for (k = 1; k < n_durations; ++k) {
        if (durations[k] < min_duration) min_duration = durations[k];
        if (durations[k] > max_duration) max_duration = durations[k];
    }
    if ((max_duration > min_period) || (min_duration < DBL_EPSILON)) return 2;

    cout << "max_duration: " << max_duration << endl ;


    // Compute the durations in terms of bin_duration
    double bin_duration = min_duration / ((double)oversample);
    int max_n_bins = (int)(ceil(max_period / bin_duration)) + oversample;

    int nthreads, blocksize = max_n_bins+1;

#pragma omp parallel
{
#if defined(_OPENMP)
    nthreads = omp_get_num_threads();
#else
    nthreads = 1;
#endif
}

    // Allocate the work arrays
    double* mean_y_0 = (double*)malloc(nthreads*blocksize*sizeof(double));
    if (mean_y_0 == NULL) {
        return -2;
    }
    double* mean_ivar_0 = (double*)malloc(nthreads*blocksize*sizeof(double));
    if (mean_ivar_0 == NULL) {
        free(mean_y_0);
        return -3;
    }


    // Pre-accumulate some factors.
    double sum_y = 0.0, sum_ivar = 0.0;
    int i;
    #pragma omp parallel for reduction(+:sum_y), reduction(+:sum_ivar)
    for (i = 0; i < N; ++i) {
        sum_y += y[i] * ivar[i];
        sum_ivar += ivar[i];
    }


    // Loop over periods and do the search
    int p;
    #pragma omp parallel for
    for (p = 0; p < n_periods; ++p) {

#if defined(_OPENMP)
        int ithread = omp_get_thread_num();
#else
        int ithread = 0;
#endif
        int block = blocksize * ithread;
        double period = periods[p];
        int n_bins = (int)(ceil(period / bin_duration)) + oversample;

        double* mean_y = mean_y_0 + block;
        double* mean_ivar = mean_ivar_0 + block;

        // This first pass bins the data into a fine-grain grid in phase from zero
        // to period and computes the weighted sum and inverse variance for each
        // bin.
        int n, ind;
        for (n = 0; n < n_bins+1; ++n) {
            mean_y[n] = 0.0;
            mean_ivar[n] = 0.0;
        }

        for (n = 0; n < N; ++n) {
            int ind = (int)(fabs(fmod(t[n], period)) / bin_duration) + 1;
            mean_y[ind] += y[n] * ivar[n];
            mean_ivar[ind] += ivar[n];
        }

        // To simplify calculations below, we wrap the binned values around and pad
        // the end of the array with the first ``oversample`` samples.
        for (n = 1, ind = n_bins - oversample; n <= oversample; ++n, ++ind) {
            mean_y[ind] = mean_y[n];
            mean_ivar[ind] = mean_ivar[n];
        }

        // To compute the estimates of the in-transit flux, we need the sum of
        // mean_y and mean_ivar over a given set of transit points. To get this
        // fast, we can compute the cumulative sum and then use differences between
        // points separated by ``duration`` bins. Here we convert the mean arrays
        // to cumulative sums.
        for (n = 1; n <= n_bins; ++n) {
            mean_y[n] += mean_y[n-1];
            mean_ivar[n] += mean_ivar[n-1];
        }

        // Then we loop over phases (in steps of n_bin) and durations and find the
        // best fit value. By looping over durations here, we get to reuse a lot of
        // the computations that we did above.
        double objective, log_like, depth, depth_err, depth_snr;
        best_objective[p] = -INFINITY;
        int k;
        for (k = 0; k < n_durations; ++k) {
            int dur = (int)(round(durations[k] / bin_duration));
            int n_max = n_bins-dur;
            for (n = 0; n <= n_max; ++n) {
                // Estimate the in-transit and out-of-transit flux
                double y_in = mean_y[n+dur] - mean_y[n];
                double ivar_in = mean_ivar[n+dur] - mean_ivar[n];
                double y_out = sum_y - y_in;
                double ivar_out = sum_ivar - ivar_in;

                // Skip this model if there are no points in transit
                if ((ivar_in < DBL_EPSILON) || (ivar_out < DBL_EPSILON)) {
                    continue;
                }

                // Normalize to compute the actual value of the flux
                y_in /= ivar_in;
                y_out /= ivar_out;

                // Either compute the log likelihood or the signal-to-noise
                // ratio
                compute_objective(y_in, y_out, ivar_in, ivar_out, obj_flag,
                        &objective, &log_like, &depth, &depth_err, &depth_snr);

                // If this is the best result seen so far, keep it
                if (y_out >= y_in && objective > best_objective[p]) {
                    best_objective[p] = objective;

                    // Compute the other parameters
                    compute_objective(y_in, y_out, ivar_in, ivar_out, (obj_flag == 0),
                            &objective, &log_like, &depth, &depth_err, &depth_snr);

                    best_depth[p]     = depth;
                    best_depth_err[p] = depth_err;
                    best_depth_snr[p] = depth_snr;
                    best_log_like[p]  = log_like;
                    best_duration[p]  = dur * bin_duration;
                    best_phase[p]     = fmod(n*bin_duration + 0.5*best_duration[p], period);
                }
            }
        }
    }

    // Clean up
    free(mean_y_0);
    free(mean_ivar_0);

    return 0;
}
*/

double weightSum(double* fluxErr, int size){
    double temp=0.0;
    for(size_t i=0; i< (size_t) size; i++){
        temp += pow((double)fluxErr[i],-2.0);
    }
    return pow(temp,-1.0);
}

double rValue(int i1, int i2, vector<double> listWeight){
    double temp=0.0;
    for(size_t i=i1; i< (size_t) i2; i++){
        temp += listWeight[i];
    }
    return temp;
}


double sValue(int i1, int i2, vector<double> listWeight, double* flux){
    double temp=0.0;
    for(size_t i=i1; i< (size_t) i2; i++){
        temp += listWeight[i]*flux[i];
    }
    return temp;
}

double dValue(vector<double> listWeight, double* flux ,double r, double s){
    double temp=0.0;
    for(size_t i=0; i<listWeight.size(); i++){
        temp += listWeight[i]* pow(flux[i],2.0);
    }
    return (double) (temp - (pow(s,2.0))) / (double) (r*( (1-r)));
}

void myBls(double* flux, double* fluxErr, double* time, int size){
    vector<double> listWeight;
    double sumW = weightSum(fluxErr, size);
    double w_i;

    for(size_t i1=0; i1<  (size_t) size; i1++){
        w_i = sumW*(pow(fluxErr[i1],-2.0));
        listWeight.push_back(w_i);
    }

    double r,s, d;

    double d_min= DBL_MAX;
    int min_i1=-1; int min_i2=-1;
    for(size_t i1=0; i1< (size_t) size; i1++){
        for(size_t i2=i1+1; i2< (size_t)size; i2++){
            r = rValue(i1,i2,listWeight);
            s = sValue(i1,i2,listWeight,flux);
            d = dValue(listWeight,flux,r,s);
            //cout << "testing i1: " << i1 << "\ti2: " << i2 << "\td: " << d  << '\n' ;
            if(d < d_min){
                d_min = d;
                min_i1=i1;
                min_i2=i2;
            }
        }
    }
   if(min_i1!=-1 && min_i2!=-1) {
       cout << "resulting i1: " << min_i1 << "\ti2: " << min_i2 << "\td: " << d_min  << '\n' ;
       cout << "Period: " << std::fixed<< time[min_i2] - time[min_i1] << '\n' ;
   }
   else cout << "could not find any pairs. latest d: " << d << endl;

}
 stapl::exit_code stapl_main(int argc, char **argv)
{
    FITS::setVerboseMode(true);
    std::string fileName = argv[1];
    try
    {

        int size = readTable(fileName);

        std::cout << "successfully readImage() \n";

        stapl::counter<stapl::default_timer> t;
        t.reset();
        t.start();
        myBls(fluxBLS, fluxErrBLS, timeBLS, size);
        double time = t.stop();
        std::cout << "finished BLS with time : " << time << "\n";
    }

    catch (FitsException&)
    // will catch all exceptions thrown by CCfits, including errors
    // found by cfitsio (status != 0)
    {
        std::cerr << " Fits Exception Thrown by test function \n";
    }

    return EXIT_SUCCESS;
}
 /*
stapl::exit_code stapl_main(int argc, char **argv)
{
    FITS::setVerboseMode(true);
    try
    {
        int size=50250;
        int n_periods=50000;
        int n_durations=5000;
        size = size - readTable("kplr008478994-2012004120508_slc.fits", size);

        std::cerr << " readImage() \n";

        //double *timeBLS, *fluxBLS, *fluxErrBLS;
        double *periods, *durations;
        periods= new double[n_periods];
        durations= new double[n_durations];

        double* best_objective = new double[n_periods];
        double* best_depth = new double[n_periods];
        double* best_depth_err = new double[n_periods];
        double* best_duration = new double[n_periods];
        double* best_phase = new double[n_periods];
        double* best_depth_snr  = new double[n_periods];
        double* best_log_like  = new double[n_periods];
        for(size_t i=0; i< (size_t) n_periods; i++){
            best_objective[i] = 0.0;
            best_depth[i] = 0.0;
            best_depth_err[i] = 0.0;
            best_duration[i] = 0.0;
            best_phase[i] = 0.0;
            best_depth_snr[i]  = 0.0;
            best_log_like[i]  = 0.0;
        }


        for(size_t i=0; i< (size_t) n_periods; i++){
            periods[i] = (rand() % 45000) + (5000);
        }
        for(size_t i=0; i< (size_t) n_durations; i++){
            durations[i] = (rand() % 5000) + (1);
        }
        sort(periods, periods + n_periods);
        sort(durations, durations + n_durations);

         stapl::counter< stapl::default_timer> t;
        t.reset();
        t.start();

        run_bls (size,timeBLS,fluxBLS,fluxErrBLS, n_periods, periods, n_durations,durations,1,1,
                 best_objective,  best_depth,  best_depth_err,  best_duration, best_phase, best_depth_snr,best_log_like);

        double read_time = t.stop();

        cout << "Results: " << std::fixed << read_time<< '\n';

        for(size_t i=0; i< (size_t) n_periods; i++){
            cout << "Curr Period: " << (int) periods[i] << endl;
            std::cout << std::fixed  << best_log_like[i] << '\t' ;
            std::cout << std::fixed  << (int) best_duration[i] << '\t';
            std::cout << std::fixed  << (int) best_phase[i] << '\n' ;
        }

    }

    catch (FitsException&)
    // will catch all exceptions thrown by CCfits, including errors
    // found by cfitsio (status != 0)
    {
        std::cerr << " Fits Exception Thrown by test function \n";
    }

    return EXIT_SUCCESS;

}
*/
