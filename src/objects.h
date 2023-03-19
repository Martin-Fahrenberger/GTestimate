#include "utils.h"
#ifndef OBJECTS_H 
#define OBJECTS_H 

class compressed_matrix {
public:
    compressed_matrix(Rcpp::RObject);
    const double* get_row(int);
    int get_ncol() const; 
    int get_nrow() const;

    bool is_row_repeated() const;
    bool is_col_repeated() const;
private:
    Rcpp::NumericMatrix mat;
    int nrow, ncol;
    bool repeat_row, repeat_col;
    std::vector<double> output;
};

class any_numeric_matrix {
public:
    any_numeric_matrix(Rcpp::RObject);
    int get_ncol() const; 
    int get_nrow() const;
    bool is_data_integer() const;

    void fill_row(int, double*);
    const Rcpp::IntegerMatrix get_raw_int() const;
    const Rcpp::NumericMatrix get_raw_dbl() const;
private:
    bool is_integer;
    int nrow, ncol;
    Rcpp::NumericMatrix dmat;
    Rcpp::IntegerMatrix imat;
};

compressed_matrix check_CM_dims(Rcpp::RObject, int, int, const char*, const char*);

Rcpp::NumericMatrix check_design_matrix(Rcpp::RObject, int);

bool check_logical_scalar(Rcpp::RObject, const char*);

int check_integer_scalar(Rcpp::RObject, const char*);

double check_numeric_scalar(Rcpp::RObject, const char*);

#endif
