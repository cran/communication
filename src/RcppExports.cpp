// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dmvnorm_cens
arma::vec dmvnorm_cens(arma::mat X, arma::rowvec mu, arma::mat Sigma, arma::uvec missingness_labels_i, std::vector< arma::uvec > nonmissing_features, bool logd, double lambda);
RcppExport SEXP _communication_dmvnorm_cens(SEXP XSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP missingness_labels_iSEXP, SEXP nonmissing_featuresSEXP, SEXP logdSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type missingness_labels_i(missingness_labels_iSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type nonmissing_features(nonmissing_featuresSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnorm_cens(X, mu, Sigma, missingness_labels_i, nonmissing_features, logd, lambda));
    return rcpp_result_gen;
END_RCPP
}
// dmvnorm_cond
arma::vec dmvnorm_cond(arma::mat X, arma::rowvec mu, arma::mat Sigma, arma::uvec labels_t, arma::uvec labels_y, bool logd, double lambda);
RcppExport SEXP _communication_dmvnorm_cond(SEXP XSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP labels_tSEXP, SEXP labels_ySEXP, SEXP logdSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type labels_t(labels_tSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type labels_y(labels_ySEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnorm_cond(X, mu, Sigma, labels_t, labels_y, logd, lambda));
    return rcpp_result_gen;
END_RCPP
}
// forward
arma::mat forward(arma::rowvec delta, arma::mat Gamma, arma::mat tstateprobs, arma::vec scale);
RcppExport SEXP _communication_forward(SEXP deltaSEXP, SEXP GammaSEXP, SEXP tstateprobsSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tstateprobs(tstateprobsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(forward(delta, Gamma, tstateprobs, scale));
    return rcpp_result_gen;
END_RCPP
}
// backward
arma::mat backward(arma::mat Gamma, arma::mat tstateprobs, arma::vec scale);
RcppExport SEXP _communication_backward(SEXP GammaSEXP, SEXP tstateprobsSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tstateprobs(tstateprobsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(backward(Gamma, tstateprobs, scale));
    return rcpp_result_gen;
END_RCPP
}
// forward_upper
arma::mat forward_upper(std::vector<arma::mat> Gammas, arma::mat tstateprobs);
RcppExport SEXP _communication_forward_upper(SEXP GammasSEXP, SEXP tstateprobsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Gammas(GammasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tstateprobs(tstateprobsSEXP);
    rcpp_result_gen = Rcpp::wrap(forward_upper(Gammas, tstateprobs));
    return rcpp_result_gen;
END_RCPP
}
// backward_upper
arma::mat backward_upper(std::vector<arma::mat> Gammas, arma::mat tstateprobs);
RcppExport SEXP _communication_backward_upper(SEXP GammasSEXP, SEXP tstateprobsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Gammas(GammasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tstateprobs(tstateprobsSEXP);
    rcpp_result_gen = Rcpp::wrap(backward_upper(Gammas, tstateprobs));
    return rcpp_result_gen;
END_RCPP
}
// hmm_cpp
Rcpp::List hmm_cpp(std::vector<arma::mat> Xs, arma::vec weights, arma::rowvec delta_init, arma::mat mus_init, std::vector<arma::mat> Sigmas_init, arma::mat Gamma_init, std::vector<arma::mat> zetas_init, std::vector< arma::uvec > nonmissing, std::vector< arma::uvec > missingness_labels, std::vector< arma::uvec > nonmissing_features, double lambda, double tol, arma::uword maxiter, double uncollapse, bool verbose, bool supervised);
RcppExport SEXP _communication_hmm_cpp(SEXP XsSEXP, SEXP weightsSEXP, SEXP delta_initSEXP, SEXP mus_initSEXP, SEXP Sigmas_initSEXP, SEXP Gamma_initSEXP, SEXP zetas_initSEXP, SEXP nonmissingSEXP, SEXP missingness_labelsSEXP, SEXP nonmissing_featuresSEXP, SEXP lambdaSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP uncollapseSEXP, SEXP verboseSEXP, SEXP supervisedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type delta_init(delta_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus_init(mus_initSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Sigmas_init(Sigmas_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma_init(Gamma_initSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type zetas_init(zetas_initSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type nonmissing(nonmissingSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type missingness_labels(missingness_labelsSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type nonmissing_features(nonmissing_featuresSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type uncollapse(uncollapseSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type supervised(supervisedSEXP);
    rcpp_result_gen = Rcpp::wrap(hmm_cpp(Xs, weights, delta_init, mus_init, Sigmas_init, Gamma_init, zetas_init, nonmissing, missingness_labels, nonmissing_features, lambda, tol, maxiter, uncollapse, verbose, supervised));
    return rcpp_result_gen;
END_RCPP
}
// hmm_autocorr_cpp
Rcpp::List hmm_autocorr_cpp(std::vector<arma::mat> Xs, arma::vec weights, arma::rowvec delta_init, arma::mat mus_init, std::vector<arma::mat> Sigmas_init, arma::mat Gamma_init, std::vector<arma::mat> zetas_init, arma::uvec labels_t, arma::uvec labels_y, double lambda, double tol, arma::uword maxiter, double uncollapse, bool verbose, bool supervised);
RcppExport SEXP _communication_hmm_autocorr_cpp(SEXP XsSEXP, SEXP weightsSEXP, SEXP delta_initSEXP, SEXP mus_initSEXP, SEXP Sigmas_initSEXP, SEXP Gamma_initSEXP, SEXP zetas_initSEXP, SEXP labels_tSEXP, SEXP labels_ySEXP, SEXP lambdaSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP uncollapseSEXP, SEXP verboseSEXP, SEXP supervisedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type delta_init(delta_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus_init(mus_initSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Sigmas_init(Sigmas_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma_init(Gamma_initSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type zetas_init(zetas_initSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type labels_t(labels_tSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type labels_y(labels_ySEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type uncollapse(uncollapseSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type supervised(supervisedSEXP);
    rcpp_result_gen = Rcpp::wrap(hmm_autocorr_cpp(Xs, weights, delta_init, mus_init, Sigmas_init, Gamma_init, zetas_init, labels_t, labels_y, lambda, tol, maxiter, uncollapse, verbose, supervised));
    return rcpp_result_gen;
END_RCPP
}
// llh_cpp
Rcpp::List llh_cpp(std::vector<arma::mat> Xs, arma::rowvec delta, arma::mat mus, std::vector<arma::mat> Sigmas_in, arma::mat Gamma, std::vector< arma::uvec > nonmissing, std::vector< arma::uvec > missingness_labels, std::vector< arma::uvec > nonmissing_features, double lambda, bool verbose);
RcppExport SEXP _communication_llh_cpp(SEXP XsSEXP, SEXP deltaSEXP, SEXP musSEXP, SEXP Sigmas_inSEXP, SEXP GammaSEXP, SEXP nonmissingSEXP, SEXP missingness_labelsSEXP, SEXP nonmissing_featuresSEXP, SEXP lambdaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Sigmas_in(Sigmas_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type nonmissing(nonmissingSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type missingness_labels(missingness_labelsSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type nonmissing_features(nonmissing_featuresSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(llh_cpp(Xs, delta, mus, Sigmas_in, Gamma, nonmissing, missingness_labels, nonmissing_features, lambda, verbose));
    return rcpp_result_gen;
END_RCPP
}
// lstateprobs_cpp
Rcpp::List lstateprobs_cpp(std::vector<arma::mat> Xs, arma::mat mus, std::vector<arma::mat> Sigmas_in, std::vector< arma::uvec > nonmissing, std::vector< arma::uvec > missingness_labels, std::vector< arma::uvec > nonmissing_features, double lambda);
RcppExport SEXP _communication_lstateprobs_cpp(SEXP XsSEXP, SEXP musSEXP, SEXP Sigmas_inSEXP, SEXP nonmissingSEXP, SEXP missingness_labelsSEXP, SEXP nonmissing_featuresSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mus(musSEXP);
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type Sigmas_in(Sigmas_inSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type nonmissing(nonmissingSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type missingness_labels(missingness_labelsSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type nonmissing_features(nonmissing_featuresSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(lstateprobs_cpp(Xs, mus, Sigmas_in, nonmissing, missingness_labels, nonmissing_features, lambda));
    return rcpp_result_gen;
END_RCPP
}
// viterbi_cpp
std::vector< std::vector<arma::uword> > viterbi_cpp(std::vector<arma::mat> lstateprobs, arma::rowvec delta, arma::mat Gamma);
RcppExport SEXP _communication_viterbi_cpp(SEXP lstateprobsSEXP, SEXP deltaSEXP, SEXP GammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<arma::mat> >::type lstateprobs(lstateprobsSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    rcpp_result_gen = Rcpp::wrap(viterbi_cpp(lstateprobs, delta, Gamma));
    return rcpp_result_gen;
END_RCPP
}
// isSilent
DataFrame isSilent(NumericVector x, double threshold);
RcppExport SEXP _communication_isSilent(SEXP xSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(isSilent(x, threshold));
    return rcpp_result_gen;
END_RCPP
}
// notSilence
DataFrame notSilence(NumericVector x, double threshold);
RcppExport SEXP _communication_notSilence(SEXP xSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(notSilence(x, threshold));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_communication_dmvnorm_cens", (DL_FUNC) &_communication_dmvnorm_cens, 7},
    {"_communication_dmvnorm_cond", (DL_FUNC) &_communication_dmvnorm_cond, 7},
    {"_communication_forward", (DL_FUNC) &_communication_forward, 4},
    {"_communication_backward", (DL_FUNC) &_communication_backward, 3},
    {"_communication_forward_upper", (DL_FUNC) &_communication_forward_upper, 2},
    {"_communication_backward_upper", (DL_FUNC) &_communication_backward_upper, 2},
    {"_communication_hmm_cpp", (DL_FUNC) &_communication_hmm_cpp, 16},
    {"_communication_hmm_autocorr_cpp", (DL_FUNC) &_communication_hmm_autocorr_cpp, 15},
    {"_communication_llh_cpp", (DL_FUNC) &_communication_llh_cpp, 10},
    {"_communication_lstateprobs_cpp", (DL_FUNC) &_communication_lstateprobs_cpp, 7},
    {"_communication_viterbi_cpp", (DL_FUNC) &_communication_viterbi_cpp, 3},
    {"_communication_isSilent", (DL_FUNC) &_communication_isSilent, 2},
    {"_communication_notSilence", (DL_FUNC) &_communication_notSilence, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_communication(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
