#ifndef UTILS_H
#define UTILS_H
//#define DEBUG

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <Rconfig.h>

#include "Rcpp.h"

#include <vector>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <algorithm>

/* Defining all R-accessible functions. */

extern "C" {

SEXP simple_good_turing (SEXP, SEXP, SEXP);

}


#endif
