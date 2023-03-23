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
#include <stapl/utility/do_once.hpp>

#include <stapl/views/array_view.hpp>
#include <stapl/containers/array/array.hpp>
#include <stapl/views/vector_view.hpp>
#include <stapl/containers/vector/vector.hpp>
#include <stapl/algorithms/algorithm.hpp>
#include <stapl/algorithms/functional.hpp>
#include <stapl/runtime.hpp>
#include <stapl/skeletons/serial.hpp>
#include <stapl/stream.hpp>
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
typedef stapl::plus<double> add_dbl_wf;

typedef std::numeric_limits< double > dbl;
const float EPSILON = 0.00000001f;
using array_type = stapl::array<double>;
using vector_type = stapl::vector<double>;

 int preprocess(stapl::vector_view<vector_type> &temp_view_flux)
{
    using domain_type = stapl::indexed_domain<size_t>;
    using vec_view_indexed_type = stapl::vector_view<stapl::vector_view<vector_type>, domain_type>;

    domain_type dom(10, 30);
    vec_view_indexed_type indFl(temp_view_flux, dom);

    cout  << "\t indFl" << endl;
    for(size_t i=indFl.domain().first() ; i< (size_t) indFl.domain().last(); i++){
        cout << '\t' <<  indFl[i] << '\n'<< flush;
    }
    cout  << "\t done " << endl<< flush;

    return 0;
}

 stapl::exit_code stapl_main(int argc, char **argv)
{
     int size = 50;
    vector_type vflux(size);
    stapl::vector_view<vector_type> view_flux(vflux);
     for(size_t i=0 ; i< (size_t) size; i++){
        view_flux[i] = 2.0;
    }
     vector_type scanned_vflux(size);
     stapl::vector_view<vector_type> scanned_view_flux(scanned_vflux);
    stapl::scan(view_flux, scanned_view_flux, add_dbl_wf(), false); // ## 4
    cout << "done prefix" << endl;
    for(size_t i=0 ; i< (size_t) size; i++){
        cout << scanned_view_flux[i] << '\t';
    }

    //int ret = preprocess(view_flux);
    //cout << temp_view_flux[12] << endl;
    //cout << ret;
    return EXIT_SUCCESS;
}