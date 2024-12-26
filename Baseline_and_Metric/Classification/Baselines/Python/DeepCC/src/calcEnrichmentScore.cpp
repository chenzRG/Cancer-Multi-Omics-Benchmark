#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double calcEnrichmentScoreCPP(IntegerVector Set, NumericVector Eso, double exponent) {
  std::vector<int> sset=Rcpp::as< std::vector<int> >(Set);
  std::vector<double> eso=Rcpp::as< std::vector<double> >(Eso);

  int N=sset.size();
  int nh=std::accumulate(sset.begin(), sset.end(), 0.0);

  if(nh == 0) {
    return 0;
  }

  double n=(-1.0)/((double)(N-nh));
  double nr=0;

  for(int j=0;j<N;j++) {
    if(sset[j]) {
      // nr += (eso[j] > 0 ? eso[j] : -eso[j]);
      nr += std::pow((eso[j] > 0 ? eso[j] : -eso[j]), exponent);
    }
  }

  double smax=0; double smin=0; double cs=0;
  for(int j=0;j<N;j++) {
    if(sset[j]) {
      // cs += (eso[j] > 0 ? eso[j] : -eso[j]) / nr;
      cs += std::pow((eso[j] > 0 ? eso[j] : -eso[j]), exponent) / nr;
    } else {
      cs += n;
    }
    if(cs>smax) {
      smax = cs;
    } else if(cs<smin) {
      smin = cs;
    }
  }
  return std::abs(smax) > std::abs(smin) ? smax : smin;
}

