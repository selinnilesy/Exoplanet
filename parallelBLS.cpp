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
//#include "benchmarks/algorithms/utilities.hpp"
#include "benchmarks/algorithms/timer.hpp"
#include <stapl/array.hpp>
#include <stapl/vector.hpp>
#include <stapl/domains/indexed.hpp>
#include <stapl/numeric.hpp>
#include <stapl/utility/do_once.hpp>
#include <stapl/algorithms/algorithm.hpp>
#include <algorithm>
#include <stapl/containers/set/set.hpp>
#include  <stapl/runtime/counter/mpi/mpi_wtime_timer.hpp>

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
typedef stapl::plus<double> double_plus;

typedef std::numeric_limits< double > dbl;
const float EPSILON = 0.00000001f;
using array_type = stapl::array<double>;
using vector_type = stapl::vector<double>;


// transform
struct powerWeight {
typedef double result_type; // ## 5
    template <typename Ref1>
    result_type operator()(Ref1 &&x) {
        cout.precision(dbl::max_digits10);
        return pow((double)x,-2.0);
  }
};
// map func
struct weightListCreate
{
private: double sumW;
public:
    weightListCreate(double s) : sumW(s) { }
    using result_type = void;
    template<typename T1, typename T2>
    result_type operator()(T1 &&val1, T2 &&output)
    {
        cout.precision(dbl::max_digits10);
        double err = val1;
        output = sumW*err;
    }
    void define_type(stapl::typer& t)
    {
        t.member(sumW);
    }
};
struct createS {
typedef double result_type; // ## 5
    template <typename Ref1, typename Ref2>
    result_type operator()(Ref1 &&flux, Ref2 &&weight ) {
        cout.precision(dbl::max_digits10);
        return flux*weight;
  }
};
struct createD
{
public:
    createD() { }
    using result_type = void;
    template<typename T1, typename T2, typename T3>
    result_type operator()(T1 &&val1, T2 &&val2, T3 &&output)
    {
        cout.precision(dbl::max_digits10);
        double flux = val1;
        double weight = val2;
        output = weight * pow(flux,2.0);
    }
};

 int preprocess(std::string fileName,  stapl::vector_view<vector_type> &fluxBLS, stapl::vector_view<vector_type> &fluxErrBLS,  stapl::vector_view<vector_type> &timeBLS)
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
    int fileSize = table.axis(1);
    cout << "read filesize: " << fileSize << endl;
    std::valarray <double> time, flux, flux_err;
    table.column("TIME").read( time, 1,fileSize );
    table.column("FLUX").read( flux, 1,fileSize );
    table.column("FLUX_ERR").read( flux_err, 1,fileSize );

    cout.precision(dbl::max_digits10);

    int toBeDeleted=0;
    for(size_t i=0; i<time.size(); i++){

        if(std::isnan(time[i]) || std::isnan(flux[i]) || std::isnan(flux_err[i]) ){
            toBeDeleted++;
        }
    }
    cout << "cleaning NaNs: " << toBeDeleted << endl;

    int newfileSize = time.size() - toBeDeleted;
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

    int numzeros=0;
    // clean 0.0s so that power does not crash the program with infs
     for(size_t i=0; i< (size_t) newfileSize; i++){
        if(((normFlux[i]<EPSILON) && (normFlux[i]>= -EPSILON)) || ((normFluxErr[i]<EPSILON) && (normFluxErr[i]>= -EPSILON)) ){
            numzeros++;
        }
    }
     cout << "cleaning also 0.0s: " << numzeros << endl;
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

void myBls( stapl::vector_view<vector_type> scannedWeights ,stapl::vector_view<vector_type> scannedWeightedFlux,
            stapl::vector_view<vector_type> time, double helper_d, int size){
    double r,s, d;
    cout.precision(dbl::max_digits10);
    int num_loc = stapl::get_num_locations();
    int loc_id = stapl::get_location_id();
    int portion = (int) ceil((1.0*size)/num_loc);
    //using domain_type = stapl::indexed_domain<size_t>;
    //using arr_view_indexed_type = stapl::array_view<stapl::array_view<array_type>, domain_type>;
    //using vec_view_indexed_type = stapl::vector_view<stapl::vector_view<vector_type>, domain_type>;
    double d_min= DBL_MAX;
    int min_i1=-1; int min_i2=-1 , upper_limit=loc_id*portion;
    cout << loc_id<< " in bls" << endl;

    if(loc_id!=num_loc-1){
        upper_limit=(loc_id+1)*portion;
    }
    else{
      upper_limit=size;
    }

    for(size_t i1=loc_id*portion; i1< (size_t) upper_limit; i1++){
        double reg1= scannedWeights[i1];
        double reg2= scannedWeightedFlux[i1];
        for(size_t i2=(i1+1); i2< (size_t) size; i2++){
            //cout << loc_id<< "-testing i1: " << i1 << "\ti2: " << i2 << "\td: " << d  << '\n' ;
            //domain_type dom(i1, i2);
            //arr_view_indexed_type v_indexedWeights(view_weight, dom);
            //r = accumulate(v_indexedWeights, 0);

            r = scannedWeights[i2] - reg1;

            //arr_view_indexed_type v_indexedS(view_weightedFlux, dom);
            //s = accumulate(v_indexedS, 0);
            s = scannedWeightedFlux[i2] - reg2;

            d = (helper_d - ( (pow(s, 2.0)) / ( 1.0 * r *  (1.0-r) )));
            d = -0.002*d + 0.032;
            if(d < d_min){
                d_min = d;
                min_i1=i1;
                min_i2=i2;
            }
        }
    }
   if(min_i1!=-1 && min_i2!=-1) {
       cout << loc_id << "'s resulting i1: " << min_i1 << "\ti2: " << min_i2 << "\td: " << d_min  << '\n' ;
       cout << loc_id << "'s Period: " << std::fixed << time[min_i2] - time[min_i1] << '\n' ;
   }
   else cout << "could not find any pairs. latest d: " << d << endl;

}
 stapl::exit_code stapl_main(int argc, char **argv)
{
    FITS::setVerboseMode(true);
    cout.precision(dbl::max_digits10);
    std::string fileName = argv[1];
    //int num_loc = stapl::get_num_locations();
    int size = 0;
    vector_type vflux;
    stapl::vector_view<vector_type> view_flux(vflux);
    vector_type vfluxerr;
    stapl::vector_view<vector_type> view_fluxerr(vfluxerr);
    vector_type vtime;
    stapl::vector_view<vector_type> view_time(vtime);

     stapl::do_once(
            [&]{
                size = preprocess(fileName, view_flux, view_fluxerr, view_time);
            });
     MPI_Bcast( &size, 1, MPI_INT, 0,  MPI_COMM_WORLD);
     int loc_id = stapl::get_location_id();


     if(loc_id==3){
         cout << "after preprocessing.. view_time \t view_flux \t view_fluxerr" << endl;
         for(size_t i=0; i< (size_t) size; i++){
            cout << std::fixed<< view_time[i] << '\t' << view_flux[i] << '\t' <<  view_fluxerr[i] << '\n';
         }
     }

    stapl::rmi_barrier();

    std::cout  << loc_id << ": successfully readImage() with size: " << size << std::endl;

    array_type squaredflux(size, 0.0);
    stapl::array_view<array_type> view_squaredfluxerr(squaredflux);
    // square errors
    stapl::transform(view_fluxerr, view_squaredfluxerr, powerWeight());
    // sum squared errors and take power of sum
    double sumW = pow(accumulate(view_squaredfluxerr, 0), -1);
    cout << "computed sumW \t" << sumW << endl << flush;

    array_type listWeight(size);
    stapl::array_view<array_type> view_weight(listWeight);
    // create weight list with sumw*err^2
    stapl::map_func(weightListCreate(sumW), view_squaredfluxerr, view_weight);

    array_type weightedFlux(size, 0.0);
    stapl::array_view<array_type> view_weightedFlux(weightedFlux);
    stapl::transform(view_flux, view_weight, view_weightedFlux, createS());

    array_type arr_d(size);
    stapl::array_view<array_type> view_d(arr_d);
    stapl::map_func(createD(), view_squaredfluxerr, view_weight, view_d);
    double helper_d = accumulate(view_d, 0);

    // prefix sum weighted flux
     vector_type scanned_weightedFlux(size);
     stapl::vector_view<vector_type> v_scanned_weightedFlux(scanned_weightedFlux);
    stapl::scan(view_weightedFlux, v_scanned_weightedFlux, double_plus(), false);
    // prefix sum weights
     vector_type scanned_weights(size);
     stapl::vector_view<vector_type> v_scanned_weights(scanned_weights);
    stapl::scan(view_weight, v_scanned_weights, double_plus(), false);

    std::vector<double> std_scanned_weights;
    std::copy(v_scanned_weights.begin(), v_scanned_weights.end(), back_inserter(std_scanned_weights));
    std::vector<double> std_scanned_weightedFlux;
    std::copy(v_scanned_weightedFlux.begin(), v_scanned_weightedFlux.end(), back_inserter(std_scanned_weightedFlux));

    /*
    using dummy = typename stapl::set<size_t>::mapper_type::dummy;

     deferred_map<block_map<id_type, location_type>,
                   decltype(map_factory), id_type, id_type>(map_factory)
     block_map<id_type, id_type>(block_size), distribution_type::blocked);
     */


    stapl::rmi_barrier();

    stapl::counter<stapl::mpi_wtime_timer> t;
    t.start();
    // do once. create time, flux and flux error.
    myBls(v_scanned_weights, v_scanned_weightedFlux, view_time, helper_d,  size);

    double time = t.stop();
    std::cout << loc_id << "-finished BLS with time : " << time << "\n";
    double max_time=0.0;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    stapl::do_once(
            [&]{ std::cout << "\nmax time : " << max_time << std::endl;}
    );

    return EXIT_SUCCESS;
}
