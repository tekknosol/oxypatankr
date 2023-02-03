#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//' Patankar scheme for modelling oxygen depletion
//'
//' This function returns a data frame containing results from the oxygen model.
//'
//' @param times An integer vector
//' @param y start concentration of O2
//' @param Flux A numeric vector of O2 fluxes. One value per iteration
//' @param Khalf A numeric vector for Khalf. One value per iteration
//' @param Theta A numeric vector for Theta. One value per iteration
//' @param Area_linear
//' @param Temp_linear
//' @param Volume_linear
//' @param iterations
//' @export
// [[Rcpp::export]]
DataFrame o2_model_lsoda_cpp(double y, NumericVector times, List forcings, List events, NumericVector Flux, NumericVector Khalf, NumericVector Theta, int iterations){
    int k = 0;
    List output_ode(1);
    NumericMatrix output(times.size(), iterations);
    NumericVector Pars;
    std::string DLLname = "o2model2";
    // mat output(times.size(), iterations);
    double it_flux, it_khalf, it_theta;
    DataFrame df;
    StringVector outnames = {"dcO2"};
    NumericVector yini = NumericVector::create(_["cO2"] = y);

    Function f("ode");
    // Function f2("dyn.load");
    // f2("inst/ext/o2model2.so");

    for(k = 0; k < iterations; k++){
        it_flux = Flux[k];
        it_khalf = Khalf[k];
        it_theta = Theta[k];

        Pars = NumericVector::create(_["Flux"] = it_flux, _["Khalf"] = it_khalf, _["Theta"] = it_theta);


        output_ode(0) = f(
        _["y"] = yini,
        _["times"] = times,
        _["func"] = "simderivs",
        _["parms"] = Pars,
        _["dllname"] = DLLname,
        _["initforc"] = "simforc",
        _["forcings"] = forcings,
        _["initfunc"] = "siminit",
        _["nout"] = 0,
        _["nspec"] = 1,
        _["outnames"] = outnames,
        _["events"] = events
    );

    NumericMatrix tmp = output_ode(0);

    output(_, k) = tmp(_,1);


    }

    return output;

}
