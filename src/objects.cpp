#include "objects.h"

compressed_matrix::compressed_matrix(Rcpp::RObject incoming) : mat(incoming) {
    if (!incoming.hasAttribute("Dims")) { 
        throw std::runtime_error("CompressedMatrix object should have 'Dims' attribute"); 
    }
    Rcpp::IntegerVector truedims=incoming.attr("Dims");
    if (truedims.size()!=2) {
        throw std::runtime_error("'Dims' attribute should be an integer vector of length 2");
    }
    nrow=truedims[0];
    ncol=truedims[1];

    if (!incoming.hasAttribute("repeat.row")) {
        throw std::runtime_error("CompressedMatrix object should have 'repeat.row' attribute");
    }
    repeat_row=check_logical_scalar(incoming.attr("repeat.row"), "'repeat.row' attribute");

    if (!incoming.hasAttribute("repeat.col")) {
        throw std::runtime_error("CompressedMatrix object should have 'repeat.col' attribute");
    }
    repeat_col=check_logical_scalar(incoming.attr("repeat.col"), "'repeat.col' attribute");
    
    // Checking dimensions.
    if (repeat_row) {
        if (mat.nrow()!=1) { 
            throw std::runtime_error("only 1 row should be repeated");
        }
    } else {
        if (mat.nrow()!=nrow) {
            throw std::runtime_error("matrix dimensions are not consistent with 'Dims'");
        }
    }

    if (repeat_col) {
        if (mat.ncol()!=1) { 
            throw std::runtime_error("only 1 col should be repeated");
        }
    } else {
        if (mat.ncol()!=ncol) {
            throw std::runtime_error("matrix dimensions are not consistent with 'Dims'");
        }
    }

    // Prefilling output vector.
    output.resize(ncol);
    if (repeat_row) {
        if (repeat_col) {
            std::fill(output.begin(), output.end(), mat[0]);
        } else {
            std::copy(mat.begin(), mat.end(), output.begin());
        }
    } 
    return;
}

const double* compressed_matrix::get_row(int index) {
    if (index>=nrow || index<0) { 
        throw std::runtime_error("requested row index out of range");
    }
    if (repeat_row) {
        return(output.data());
    } 
    if (repeat_col) {
        std::fill(output.begin(), output.end(), *(mat.begin()+index));
    } else { 
        auto cmat=mat.row(index);
        std::copy(cmat.begin(), cmat.end(), output.begin());
    }
    return(output.data());
}

bool compressed_matrix::is_row_repeated () const {
    return repeat_row;
}

bool compressed_matrix::is_col_repeated () const {
    return repeat_col;
}

int compressed_matrix::get_ncol() const { return ncol; }

int compressed_matrix::get_nrow() const { return nrow; }

/* Methods for any numeric matrix */

any_numeric_matrix::any_numeric_matrix(Rcpp::RObject incoming) : is_integer(incoming.sexp_type()==INTSXP) {
    if (is_integer){ 
        imat=Rcpp::IntegerMatrix(incoming);
        nrow=imat.nrow();
        ncol=imat.ncol();
    } else {
        dmat=Rcpp::NumericMatrix(incoming);
        nrow=dmat.nrow();
        ncol=dmat.ncol();
    }
    return;
}

void any_numeric_matrix::fill_row(int index, double* ptr) {
    if (is_integer) {
        auto current=imat.row(index);
        std::copy(current.begin(), current.end(), ptr);
    } else {
        auto current=dmat.row(index);
        std::copy(current.begin(), current.end(), ptr);;
    }
    return;
}

bool any_numeric_matrix::is_data_integer () const {
    return is_integer;
}

const Rcpp::IntegerMatrix any_numeric_matrix::get_raw_int() const {
    return imat;
}

const Rcpp::NumericMatrix any_numeric_matrix::get_raw_dbl() const {
    return dmat;
}

int any_numeric_matrix::get_ncol() const { return ncol; }

int any_numeric_matrix::get_nrow() const { return nrow; }

/* Methods to check the dimensions of any object. */

compressed_matrix check_CM_dims(Rcpp::RObject incoming, int nrow, int ncol, const char* current, const char* ref) {
     compressed_matrix out(incoming);
     if (out.get_nrow()!=nrow || out.get_ncol()!=ncol) {
         std::stringstream err;
         err << current << " and " << ref << " matrices do not have the same dimensions";
         throw std::runtime_error(err.str().c_str());
     }
     return out;
}

Rcpp::NumericMatrix check_design_matrix(Rcpp::RObject design, int nlibs) {
    Rcpp::NumericMatrix X(design);
    if (X.nrow()!=nlibs) {
        throw std::runtime_error("number of rows in the design matrix should be equal to the number of libraries");
    }
    return X;
}

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

bool check_logical_scalar(Rcpp::RObject x, const char* thing) {
    return check_scalar_value<bool, Rcpp::LogicalVector>(x, "logical", thing);
}

int check_integer_scalar(Rcpp::RObject x, const char* thing) {
    return check_scalar_value<int, Rcpp::IntegerVector>(x, "integer", thing);
}

double check_numeric_scalar(Rcpp::RObject x, const char* thing) {
    return check_scalar_value<double, Rcpp::NumericVector>(x, "double-precision", thing);
}

