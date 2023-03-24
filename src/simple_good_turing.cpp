#include "utils.h"
#include "Rcpp.h"

//' Implements the simple version of the Good-Turing frequency estimator in C++.
//' This is based on the C code written by Geoffrey Sampson from Sussex University.
//' It takes in a vector of observed frequencies and another vector of the same
//' length of frequencies (of observed frequencies). The first vector must be
//' sorted in ascending order. It also takes a numeric scalar which describes
//' the confidence factor.
//'
//' @param obs: observed values (must be sorted in ascending order)
//' @param freq: frequencies of observed values
//' @param conf: confidence factor for internal fit
//' @return named list

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

SEXP simple_good_turing (SEXP obs, SEXP freq, SEXP conf) {
    BEGIN_RCPP

    Rcpp::IntegerVector Obs(obs);
    Rcpp::IntegerVector Freq(freq);
    const int nrows=Obs.size();
    if (nrows!=Freq.size()) {
        throw std::runtime_error("lengths of obs and freq vectors must match");
    }

    // Checking the confidence factor.
    double confid_factor=check_numeric_scalar(conf, "confidence factor");

	// Computing constant values.
	double bigN=0;
	double XYs=0, meanX=0, meanY=0, Xsquares=0;
    std::vector<double> log_obs(nrows);
	const int last=nrows-1;

    for (int i=0; i<nrows; ++i) {
		const int& o=Obs[i];
        const int& f=Freq[i];
		bigN+=o*f;

		// Computing log data.
		const int& x=(i==0 ? 0 : Obs[i-1]);
		const double& logO=(log_obs[i]=std::log(double(o)));
		const double logZ=std::log(2*f/double(i==last ? 2*(o-x) : Obs[i+1]-x));
		meanX+=logO;
		meanY+=logZ;
		XYs+=logO*logZ;
		Xsquares+=logO*logO;
	}

	meanX/=nrows;
	meanY/=nrows;
	XYs-=meanX*meanY*nrows;
	Xsquares-=meanX*meanX*nrows;
	const double slope=XYs/Xsquares;

	// Computing other various bits and pieces.
	const double& PZero = ((nrows==0 || Obs[0]!=1) ? 0 : Freq[0]/double(bigN));

    // Collecting results.
    double bigNprime=0;
	bool indiffValsSeen=false;
    Rcpp::NumericVector outp(nrows);

    for (int i=0; i<nrows; ++i) {
        const int& o=Obs[i];
        const int& f=Freq[i];
        double& op=outp[i];

		const int next_obs=o+1;
		const double y = next_obs*std::exp(slope*(std::log(double(next_obs))-log_obs[i])); // don't need intercept, cancels out.
		if (i==last || Obs[i+1]!=next_obs) {
            indiffValsSeen=true;
        }

		if (!indiffValsSeen) {
            const int& next_n=Freq[i+1];
            const double x = next_obs*next_n/double(f);
            if (std::abs(x - y) <= confid_factor * x * std::sqrt(1.0/next_n + 1.0/double(f))) { // Simplified expression.
                indiffValsSeen=true;
            } else {
                op=x;
			}
		}
		if (indiffValsSeen) {
			op=y;
        }
        bigNprime+=op*f;
	}

	// Running through them to compute the remaining bit.
	const double factor=(1.0-PZero)/bigNprime;
    for (auto& op : outp) {
        op*=factor;
    }

    //return Rcpp::List::create(Rcpp::NumericVector::create(PZero), outp);
    return outp;
    END_RCPP
}


