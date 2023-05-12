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
#include <stapl/numeric.hpp>
#include <stapl/utility/do_once.hpp>
#include <stapl/algorithms/algorithm.hpp>
#include <algorithm>
#include <stapl/containers/distribution/distribution.hpp>
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
using array_type = stapl::array< double>;
using vector_type = stapl::vector<double>;
using vector_type_int = stapl::vector<int>;
using nested_array_type = stapl::array<vector_type>;
using nested_array_type_int = stapl::array< std::vector<int>>;

#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif
// transform
struct powerWeight {
typedef double result_type; // ## 5
    template <typename Ref1>
    result_type operator()(Ref1 &&x) {
        cout.precision(dbl::max_digits10);
        return 1.0/(x*x);
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
        double sf = val1;
        double weight = val2;
        output = weight * (1.0 / sf);
    }
};
struct BLS_mapfunc
{
public:
    int l_size;
    int prev_size;
    int res_i1;
    int res_i2;
    double res_mind;
    double helper_d;
    BLS_mapfunc(int l_size, int prev_size, double helper_d) :  l_size(l_size), prev_size(prev_size), res_i1(-1), res_i2(-1), res_mind(DBL_MAX), helper_d(helper_d)
    {}
    using result_type = void;
    template<typename T1, typename T2>
    result_type operator()(T1 &&scannedWeights, T1 &&scannedWeightedFlux, T2 &&sizes)
    {
        double r,s, d;
        cout.precision(dbl::max_digits10);
        int num_loc = stapl::get_num_locations();
        int loc_id = stapl::get_location_id();
        int next_size=0;

        if(loc_id != num_loc-1)
        {
            auto future_weights = scannedWeights.container().get_element_split(loc_id+1);
            auto future_fluxedWeights = scannedWeightedFlux.container().get_element_split(loc_id+1);
            next_size= sizes[loc_id+1];

            for(size_t i1=0; i1< (size_t) l_size; i1++){
                double reg1= scannedWeights[i1];
                double reg2= scannedWeightedFlux[i1];
                for(size_t i2=(i1+1); i2< (size_t) l_size; i2++){
                    r = scannedWeights[i2] - reg1;
                    s = scannedWeightedFlux[i2] - reg2;
                    d = (helper_d - ( s*s)) / ( 1.0 * r *  (1.0-r) );
                    if(d < res_mind){
                        res_mind  = d;
                        res_i1= prev_size + i1;
                        res_i2= prev_size + i2;
                    }
                }
            }

            std::vector<double> next_weights = future_weights.get();
            std::vector<double> next_fluxedWeights = future_fluxedWeights.get();

            for(int i=loc_id+2; i<num_loc; i++){

                auto future_weights = scannedWeights.container().get_element_split(i);
                auto future_fluxedWeights = scannedWeightedFlux.container().get_element_split(i);

                for(size_t i1=0; i1< (size_t) l_size; i1++){
                    double reg1= next_weights[i1];
                    double reg2= next_fluxedWeights[i1];
                    for(size_t i2= 0; i2< (size_t) next_size; i2++){
                        r = next_weights[i2] - reg1;
                        s = next_fluxedWeights[i2] - reg2;
                        d = (helper_d - ( s*s)) / ( 1.0 * r *  (1.0-r) );
                        if(d < res_mind){
                            res_mind  = d;
                            res_i1= prev_size + i1;
                            res_i2= prev_size + i2;
                        }
                    }
                }
                next_weights = future_weights.get();
                next_fluxedWeights = future_fluxedWeights.get();
                next_size =  sizes[i];
            }
            for(size_t i1=0; i1< (size_t) l_size; i1++){
                double reg1= next_weights[i1];
                double reg2= next_fluxedWeights[i1];
                for(size_t i2=0; i2< (size_t) next_size; i2++){
                    r = next_weights[i2] - reg1;
                    s = next_fluxedWeights[i2] - reg2;
                    d = (helper_d - ( s*s)) / ( 1.0 * r *  (1.0-r) );
                    if(d < res_mind){
                        res_mind  = d;
                        res_i1= prev_size + i1;
                        res_i2= prev_size + i2;
                    }
                }
            }
        }
    }

    void define_type(stapl::typer& t)
    {
        t.member(l_size);
        t.member(prev_size);
        t.member(res_i1);
        t.member(res_i2);
        t.member(res_mind);
        t.member(helper_d);
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

 typedef struct {
    double d;
    int min_i1;
    int min_i2;
 } solution_t;

 /* custom structs and reduction, taken from here
  * https://stackoverflow.com/questions/65874060/an-efficient-way-to-perform-an-all-reduction-in-mpi-of-a-value-based-on-another
  */
 void compute_min_struct(void *in, void *inout, int *len, MPI_Datatype *type){
    solution_t *invals    = (solution_t *) in;
    solution_t *inoutvals = (solution_t *) inout;
    for (int i=0; i<(*len); i++){
        if(inoutvals[i].d < invals[i].d){
            inoutvals[i].d  = invals[i].d;
            inoutvals[i].min_i1  = invals[i].min_i1;
            inoutvals[i].min_i2  = invals[i].min_i2;
        }
    }
}
void defineStruct(MPI_Datatype *tstype) {
    const int count = 3;
    int          blocklens[count];
    MPI_Datatype types[count];
    MPI_Aint     disps[count];

    for (int i=0; i < count; i++) {
        if(i==0){
            types[i] = MPI_DOUBLE;
        }
        else{
            types[i] = MPI_INT;
        }
        blocklens[i] = 1;
    }
    disps[0] = offsetof(solution_t,d);
    disps[1] = offsetof(solution_t,min_i1);
    disps[2] = offsetof(solution_t,min_i2);

    MPI_Type_create_struct(count, blocklens, disps, types, tstype);
    MPI_Type_commit(tstype);
}

void myBls( stapl::vector_view<vector_type> scannedWeights ,stapl::vector_view<vector_type> scannedWeightedFlux,
            stapl::vector_view<vector_type> time, double helper_d, size_t size){
    double r,s, d;
    cout.precision(dbl::max_digits10);
    int num_loc = stapl::get_num_locations();
    int loc_id = stapl::get_location_id();

    MPI_Datatype structtype;
    MPI_Op       minop_struct;
    defineStruct(&structtype);
    MPI_Op_create(compute_min_struct, 1, &minop_struct);

    solution_t local_solution = { DBL_MAX, -1,-1 };
    solution_t global_solution = { DBL_MAX, -1,-1 };
    double max_time=0.0;

    stapl::counter<stapl::mpi_wtime_timer> t;
    t.start();

    size_t grid= size/sqrt(num_loc);
    size_t start = size - sqrt(num_loc - loc_id) * grid;
    size_t end = size - sqrt(num_loc - loc_id -1) * grid;


    for(size_t i1=0; i1< (size_t) end-start; i1++){
        double reg1= scannedWeights[i1];
        double reg2= scannedWeightedFlux[i1];

        for(size_t i2=(i1+1); i2< (size_t) end-start; i2++){
            r = scannedWeights[i2] - reg1;
            s = scannedWeightedFlux[i2] - reg2;
            d = (helper_d - ( s*s)) / ( 1.0 * r *  (1.0-r) );
            if(d < local_solution.d){
                local_solution.d  = d;
                local_solution.min_i1=i1;
                local_solution.min_i2=i2;
            }
        }
    }

    MPI_Reduce(&local_solution, &global_solution, 1, structtype, minop_struct, 0, MPI_COMM_WORLD);
    double measured_time = t.stop();
    MPI_Reduce(&measured_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    stapl::do_once(
    [&]{
        std::cout << "\nmax time : " << max_time << std::endl;
        if(global_solution.min_i1!= -1 && global_solution.min_i2!= -1) {
            cout << " Resulting i1: " << global_solution.min_i1 << "\ti2: " << global_solution.min_i2 << "\td: " << global_solution.d  << '\n' ;
            cout <<" Period: " << std::fixed << time[global_solution.min_i2] - time[global_solution.min_i1] << '\n' ;
        }
        else cout << "could not find any pairs. latest d: " << d << endl;
    });
}
 stapl::exit_code stapl_main(int argc, char **argv)
{
    FITS::setVerboseMode(true);
    cout.precision(dbl::max_digits10);
    std::string fileName = argv[1];
    int num_loc = stapl::get_num_locations();
    size_t size = 0;
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

    stapl::rmi_barrier();

    std::cout  << loc_id << ": successfully readImage() with size: " << size << std::endl;

    array_type squaredflux((size), 0.0);
    stapl::array_view<array_type> view_squaredfluxerr(squaredflux);
    // square errors
    stapl::transform(view_fluxerr, view_squaredfluxerr, powerWeight());
    // sum squared errors and take power of sum
    double sumW = pow(accumulate(view_squaredfluxerr, 0), -1);
    cout << "computed sumW \t" << sumW << endl << flush;

    array_type listWeight((size));
    stapl::array_view<array_type> view_weight(listWeight);
    // create weight list with sumw*err^2
    stapl::map_func(weightListCreate(sumW), view_squaredfluxerr, view_weight);

    array_type weightedFlux((size), 0.0);
    stapl::array_view<array_type> view_weightedFlux(weightedFlux);
    stapl::transform(view_flux, view_weight, view_weightedFlux, createS());

    array_type arr_d((size));
    stapl::array_view<array_type> view_d(arr_d);
    stapl::map_func(createD(), view_squaredfluxerr, view_weight, view_d);
    double helper_d = accumulate(view_d, 0);

    // prefix sum weighted flux
     vector_type scanned_weightedFlux((size));
     stapl::vector_view<vector_type> v_scanned_weightedFlux(scanned_weightedFlux);
    stapl::scan(view_weightedFlux, v_scanned_weightedFlux, double_plus(), false);
    // prefix sum weights
     vector_type scanned_weights((size));
     stapl::vector_view<vector_type> v_scanned_weights(scanned_weights);
    stapl::scan(view_weight, v_scanned_weights, double_plus(), false);

    /*
    std::vector<double> std_scanned_weights;
    std::copy(v_scanned_weights.begin(), v_scanned_weights.end(), back_inserter(std_scanned_weights));
    std::vector<double> std_scanned_weightedFlux;
    std::copy(v_scanned_weightedFlux.begin(), v_scanned_weightedFlux.end(), back_inserter(std_scanned_weightedFlux));
     */
    /*
    using dummy = typename stapl::set<size_t>::mapper_type::dummy;

     deferred_map<block_map<id_type, location_type>,
                   decltype(map_factory), id_type, id_type>(map_factory)
     block_map<id_type, id_type>(block_size), distribution_type::blocked);
     */

     /*
     * vector_type m_scanned_weights;
     stapl::vector_view<vector_type> mv_scanned_weights(m_scanned_weights);
     vector_type m_scanned_weightedFlux;
     stapl::vector_view<vector_type> mv_scanned_weightedFlux(m_scanned_weightedFlux);
     vector_type m_vtime;
     stapl::vector_view<vector_type> mv_time(m_vtime);
    for(size_t i=start; i< end; i++){
        mv_scanned_weights.add(v_scanned_weights[i]);
        mv_scanned_weightedFlux.add(v_scanned_weightedFlux[i]);
        mv_time.add(vtime[i]);
    }
    mv_scanned_weights.container().distribution().synchronize_metadata();
    mv_scanned_weightedFlux.container().distribution().synchronize_metadata();
    mv_time.container().distribution().synchronize_metadata();
     */


    size_t grid= size/sqrt(num_loc);
    size_t start = size - sqrt(num_loc - loc_id) * grid;
    size_t end = size - sqrt(num_loc - loc_id -1) * grid;

    //std::vector<double> w, wf;
    vector_type w;
    stapl::vector_view<vector_type> v_w(w);
    vector_type wf;
    stapl::vector_view<vector_type> v_wf(wf);
    for(size_t i=start; i< end; i++){
        v_w.add(v_scanned_weights[i]);
        v_wf.add(v_scanned_weightedFlux[i]);
    }
    nested_array_type w_location_vecs(num_loc);
    stapl::array_view<nested_array_type> v_w_location_vecs(w_location_vecs);
    v_w_location_vecs[loc_id] = w;
    nested_array_type wf_location_vecs(num_loc);
    stapl::array_view<nested_array_type> v_wf_location_vecs(wf_location_vecs);
    v_wf_location_vecs[loc_id] = wf;

    std::vector<int> size_locations(num_loc);
    size_locations[loc_id] = end-start;

    nested_array_type_int size_location_vecs(num_loc);
    stapl::array_view<nested_array_type_int> v_size_location_vecs(size_location_vecs);
    v_size_location_vecs[loc_id] = size_locations;

    int lsize = end-start;
    int prevsize=0;
    MPI_Exscan(&lsize, &prevsize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    //cout << loc_id<< "'s read: mv_scanned_weights-" << mv_scanned_weights[5] << "\tmv_scanned_weightedFlux-" << mv_scanned_weightedFlux[10]  << '\n' ;

    cout << loc_id << ": came to barrier \t" << endl << flush;
    stapl::rmi_barrier();
    cout << loc_id << ": starting myBLS \t" << endl << flush;
    //myBls(mv_scanned_weights, mv_scanned_weightedFlux, mv_time, helper_d,  size);
    stapl::map_func(BLS_mapfunc(end-start,prevsize, helper_d), v_w_location_vecs, v_wf_location_vecs, v_size_location_vecs);

    return EXIT_SUCCESS;
}
