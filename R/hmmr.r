#' Viterbi algorithm implemented in R
#'
#' This function implements a Viterbi algorithm, given observed outcomes, transition and emission matrices, and initial probabilities.
#'
#' @param obs size T integer vector of observed outcomes, values from 1 to m
#' @param transMtx n x n transition matrix from row to column
#' @param emisMtx n x m emission matrix of observed data (column) given states (row)
#' @param pi size n vector of initial state probabilities
#'
#' @return a list containing maximum likelihood (ml) and the Viterbi path (path)
#'
#' @examples
#' obs <- c(1,3,2) 
#' A <- matrix(c(0.8,0.2,0.4,0.6),2,2,byrow=TRUE)
#' B <- matrix(c(0.88,0.10,0.02,0.1,0.6,0.3),2,3,byrow=TRUE)
#' pi <- c(0.7,0.3)
#' viterbiHMMr(obs, A, B, pi)
#' 
#' @export
viterbiHMMr <- function(obs, transMtx, emisMtx, pi) {
    T <- length(obs)
    ns <- nrow(transMtx)
    delta <- matrix(NA,ns,T)
    phi <- matrix(NA,ns,T)
    delta[,1] <- pi * emisMtx[,obs[1]]
    for(t in seq(2,T,1)) {
        for(i in 1:ns) { ## no need to use double loop
            v <- delta[,t-1] * transMtx[,i]
            delta[i,t] <- max(v) * emisMtx[i,obs[t]]
            phi[i,t] <- which.max(v)
        }
    }
    path <- vector(length=T)
    ml <- max(delta[,T])    
    path[T] <- which.max(delta[,T])
    for(t in seq(T-1,1,-1)) {
        path[t] = phi[path[t+1],t+1]
    }
    return(list(ml=ml,path=path,delta=delta,phi=phi))
}

#' Forward-backward algorithm implemented in R
#'
#' This function implements a forward-backward algorithm for a generic HMM in R, given observed outcomes, transition and emission matrices, and initial probabilities.
#'
#' @param obs size T integer vector of observed outcomes, values from 1 to m
#' @param transMtx n x n transition matrix from row to column
#' @param emisMtx n x m emission matrix of observed data (column) given states (row)
#' @param pi size n vector of initial state probabilities
#'
#' @return a matrix containing conditional probability of each possible states given the observed outcomes
#'
#' @examples
#' obs <- c(1,3,2) 
#' A <- matrix(c(0.8,0.2,0.4,0.6),2,2,byrow=TRUE)
#' B <- matrix(c(0.88,0.10,0.02,0.1,0.6,0.3),2,3,byrow=TRUE)
#' pi <- c(0.7,0.3)
#' forwardBackwardHMMr(obs, A, B, pi)
#' 
#' @export
forwardBackwardHMMr <- function(obs, transMtx, emisMtx, pi) {
    T <- length(obs)
    ns <- nrow(transMtx)
    alpha <- matrix(0, ns, T)
    beta <- matrix(0, ns, T)
    alpha[,1] <- pi * emisMtx[,obs[1]]
    for(t in seq(2,T,1)) {
        for(i in 1:ns) {
            alpha[i,t] = sum(alpha[,t-1] * transMtx[,i]) * emisMtx[i,obs[t]]
        }
    }
    beta[,T] <- 1
    for(t in seq(T-1,1,-1)) {
        for(i in 1:ns) {
            beta[i,t] = sum(beta[,t+1] * transMtx[i,] * emisMtx[,obs[t+1]])
        }
    }
    condProb <- alpha * beta
    condProb <- condProb / matrix(colSums(condProb),ns,T,byrow=TRUE)
    return(condProb)
}
