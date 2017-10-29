#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

void forwardLoop(IntegerVector& obs, NumericMatrix& transMtx,
               NumericMatrix& emisMtx, NumericVector& pi,
               NumericMatrix& alpha) {
  int T = (int)obs.size();  
  int ns = (int)pi.size();
  
  for(int i=0; i < ns; ++i) 
    alpha(i,0) = pi(i) * emisMtx(i, obs[0]-1);
  
  for(int t=1; t < T; ++t) {
    for(int i=0; i < ns; ++i) {
      alpha(i,t) = 0;
      for(int j=0; j < ns; ++j) 
        alpha(i,t) += ( alpha(j,t-1) * transMtx(j,i) );
      alpha(i,t) *= emisMtx(i, obs[t]-1); // why did I do this?
    }
  }
}

void backwardLoop(IntegerVector& obs, NumericMatrix& transMtx,
               NumericMatrix& emisMtx, NumericVector& pi,
               NumericMatrix& beta) {
  int T = (int)obs.size();  
  int ns = (int)pi.size();
  
  for(int i=0; i < ns; ++i) 
    beta(i,T-1) = 1;
  
  for(int t=T-2; t >=0; --t) {
    for(int i=0; i < ns; ++i) {
      beta(i,t) = 0;
      for(int j=0; j < ns; ++j) 
        beta(i,t) += ( beta(j,t+1) * transMtx(i,j) * emisMtx(j, obs[t+1]-1) );
    }
  }
}

//' Forward-backward algorithm implemented in Rcpp
//'
//' This function implements a forward-backward algorithm of a generic HMM, given observed outcomes, transition and emission matrices, and initial probabilities.
//'
//' @param obs size T integer vector of observed outcomes, values from 1 to m
//' @param transMtx n x n transition matrix from row to column
//' @param emisMtx n x m emission matrix of observed data (column) given states (row)
//' @param pi size n vector of initial state probabilities
//' @return matrix of conditional probability of each state given the observed outcomes
//' @examples 
//' obs <- c(1,3,2)
//' A <- matrix(c(0.8,0.2,0.4,0.6),2,2,byrow=TRUE)
//' B <- matrix(c(0.88,0.10,0.02,0.1,0.6,0.3),2,3,byrow=TRUE)
//' pi <- c(0.7,0.3)
//' forwardBackwardHMMc(obs, A, B, pi)
//' @export
// [[Rcpp::export]]
NumericMatrix forwardBackwardHMMc(IntegerVector obs, NumericMatrix transMtx, NumericMatrix emisMtx, NumericVector pi) {
  int T = (int)obs.size();  
  int ns = (int)pi.size();
  NumericMatrix alpha(ns, T);
  NumericMatrix beta(ns, T);
  
  forwardLoop(obs, transMtx, emisMtx, pi, alpha);
  backwardLoop(obs, transMtx, emisMtx, pi, beta);
  NumericMatrix condProb(ns, T);
  
  for(int t=0; t < T; ++t) {
    double sum = 0;
    
    for(int i=0; i < ns; ++i) 
      sum += ( alpha(i,t) * beta(i,t) );
    for(int i=0; i < ns; ++i) 
      condProb(i,t) = alpha(i,t) * beta(i,t) / sum;
  }
  return condProb;
}

double viterbiLoop(IntegerVector& obs, NumericMatrix& transMtx,
                   NumericMatrix& emisMtx, NumericVector& pi,
                   NumericMatrix& delta, IntegerMatrix& phi) {
  int T = (int)obs.size();  
  int ns = (int)pi.size();
  for(int i=0; i < ns; ++i) 
    delta(i,0) = pi(i) * emisMtx(i, obs[0]-1);
  for(int t=1; t < T; ++t) {
    for(int i=0; i < ns; ++i) {
      for(int j=0; j < ns; ++j) {
        double v = delta(j,t-1) * transMtx(j, i);
        if ( v > delta(i,t) ) {
          delta(i,t) = v;
          phi(i,t) = j;
        }
      }
      delta(i,t) *= emisMtx(i, obs[t]-1);
    }
  }
  double ml = -1;
  for(int i=0; i < ns; ++i) {
    if ( delta(i,T-1) > ml ) ml = delta(i,T-1);
  }
  return ml;
}

//' Viterbi algorithm implemented in Rcpp
//'
//' This function implements a Vitervi algorithm in generic HMM, given observed outcomes, transition and emission matrices, and initial probabilities.
//'
//' @param obs size T integer vector of observed outcomes, values from 1 to m
//' @param transMtx n x n transition matrix from row to column
//' @param emisMtx n x m emission matrix of observed data (column) given states (row)
//' @param pi size n vector of initial state probabilities
//' @return a matrix containing conditional probability of each possible states given the observed outcomes
//' @examples 
//' obs <- c(1,3,2)
//' A <- matrix(c(0.8,0.2,0.4,0.6),2,2,byrow=TRUE)
//' B <- matrix(c(0.88,0.10,0.02,0.1,0.6,0.3),2,3,byrow=TRUE)
//' pi <- c(0.7,0.3)
//' viterbiHMMc(obs, A, B, pi)
//' @export
// [[Rcpp::export]]
List viterbiHMMc(IntegerVector obs, NumericMatrix transMtx, NumericMatrix emisMtx, NumericVector pi) {
  int T = (int)obs.size();  
  int ns = (int)pi.size();
  if ( ( ns != transMtx.ncol() ) || ( ns != transMtx.nrow() ) )
    stop("transMtx must be square and have same dimension with pi");
  if ( ns != emisMtx.nrow() )
    stop("emisMtx must have the same number of rows to transMtx");

  NumericMatrix delta(ns, T);
  IntegerMatrix phi(ns, T);
  std::fill(delta.begin(), delta.end(), -1);
  
  viterbiLoop(obs, transMtx, emisMtx, pi, delta, phi);
  double ml = -1;
  IntegerVector paths(T);
  for(int i=0; i < ns; ++i) {
    if ( ml < delta(i,T-1) ) {
      ml = delta(i,T-1);
      paths[T-1] = i;
    }
  }  
  for(int i=T-1; i > 0; --i) {
    paths[i-1] = phi(paths[i],i);
  }
  
  return ( List::create(Named("ml") = ml,
                        Named("path") = paths,
                        Named("delta") = delta,
                        Named("phi") = phi
  ) );
}
