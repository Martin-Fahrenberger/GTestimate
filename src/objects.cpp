#include "objects.h"

template<typename T, class V>
T check_scalar_value (Rcpp::RObject val, const char* type, const char* thing) {
    V x(val);
    if (x.size()!=1) {
        std::stringstream err;
        err << "expected " << type << " scalar for the " << thing;
        throw std::runtime_error(err.str().c_str());
    }
    return x[0];
}

double check_numeric_scalar(Rcpp::RObject x, const char* thing) {
    return check_scalar_value<double, Rcpp::NumericVector>(x, "double-precision", thing);
}

