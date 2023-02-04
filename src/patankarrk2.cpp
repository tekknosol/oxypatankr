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
DataFrame o2_model_patankarrk2_cpp(NumericVector times, NumericVector y, NumericVector onset, NumericVector Flux, NumericVector Khalf, NumericVector Theta, NumericVector Area_linear, NumericVector Temp_linear, NumericVector Volume_linear, int iterations) {
  int n_times = times.size();
  mat dcO2(n_times, iterations);
  rowvec cO20(iterations, fill::value(y(0)));
  //int len_y0 =  y.size();
  vec r(1);
  int i = 0, k = 0;
  int dt = 1;
  vec MichaelisMenten(1), ArrheniusCorrection(1), cO2_vec(1);
  vec SedimentFlux(1);
  // double p0 = 1e-10, d0, ydat, c_rep, p1 = 1e-10, d1, p = 0.5 * (p0 + p1), d = 0.5 * (p0 + p1), cproxy_d;
  double p0 = 1e-10, d0, ydat, c_rep, p1 = 1e-10, d1, p, d, cproxy_d;
  double it_flux, it_khalf, it_theta;
  vec cproxy;
  // mat eye(len_y0, len_y0, fill::ones);
  // mat eyetilde(len_y0, len_y0, fill::zeros);
  mat avec(1, 1);
  vec mean_cO2(n_times), median_cO2(n_times), sd_cO2(n_times), up_cO2(n_times), lp_cO2(n_times);
  vec P;

  DataFrame df;

  dcO2.row(0) = cO20;

  p = 0.5 * (p0 + p1);
  d = 0.5 * (p0 + p1);


  for (k = 0; k < iterations; k++){
    it_flux = Flux[k];
    it_khalf = Khalf[k];
    it_theta = Theta[k];


    for (i = 1; i < n_times; i++){
      if (onset[i] == 1)
      {
        dcO2(i,k) = y[i];

      } else
      {

        SedimentFlux[0] = Area_linear[i] * it_flux ;
        MichaelisMenten[0] = ((dcO2(i - 1, k)) / (it_khalf + dcO2(i - 1, k)));
        ArrheniusCorrection[0] = pow(it_theta, (Temp_linear[i] - 20));

         p0 = 1e-10;
        d0 = std::abs(SedimentFlux[0] * MichaelisMenten[0] * ArrheniusCorrection[0] / Volume_linear[i]);


        ydat = dcO2(i - 1, k);

        avec.fill(dt * d0 / ydat + 1);

        c_rep = ydat;


        r(0) = ydat + dt * p0;

        cproxy = solve(avec, r);

        cproxy_d = cproxy(0,0);

        SedimentFlux[0] = Area_linear[i] * it_flux;
        MichaelisMenten[0] = ((cproxy_d) / (it_khalf + cproxy_d));
        ArrheniusCorrection[0] = pow(it_theta, (Temp_linear[i] - 20));

         p1 = 1e-10;
        d1 =  std::abs(SedimentFlux[0] * MichaelisMenten[0] * ArrheniusCorrection[0] / Volume_linear[i]);

         p = 0.5 * (p0 + p1);
         d = 0.5 * (d0 + d1);

        avec.fill(dt * d / cproxy_d + 1);

        r(0) = ydat + dt * p;

        cO2_vec = solve(avec, r);


        dcO2(i,k) = cO2_vec(0,0);


      }


    }

  }

  mean_cO2 = mean(dcO2, 1);
  median_cO2 = median(dcO2, 1);
  sd_cO2 = stddev(dcO2, 0, 1);
  P = {0.975};
  up_cO2 = quantile(dcO2, P, 1);
  P = {0.025};
  lp_cO2 = quantile(dcO2, P, 1);

  df = DataFrame::create(
    Named("oxygen_mean") = mean_cO2,
    Named("oxygen_median") = median_cO2,
    Named("oxygen_sd") = sd_cO2,
    Named("oxygen_upperPercentile") = up_cO2,
    Named("oxygen_lowerPercentile") = lp_cO2
  );
  return(df);


}
