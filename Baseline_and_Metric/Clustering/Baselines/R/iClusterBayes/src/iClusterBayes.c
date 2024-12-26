/* Qianxing Mo (qianxing.mo@moffitt.org), Department of Biostatistics & Bioinformatics, H. Lee Moffitt Cancer Center and Research Institute */
/* Code programs for iClusterBayes
   last updated 11/05/2015, change logp to logpnull, logq to logqnull
*/
#define USE_FC_LEN_T
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h> 
#include <R_ext/Lapack.h>
#include <Rconfig.h>
/* Note include "iClusterPlus.h" must be put behind include <math.h> */
/* Otherwise beta variable will confict with #define beta Rf_beta in <math.h> */
#include "iClusterPlus.h"

#ifndef FCONE
#define FCONE
#endif

/* x = x + y */
void ivadd(int * x, int *y, int size) {
  int i;
  for(i = 0; i < size; i++) {
    x[i] = x[i] + y[i];
  }
}
/*scale(Z,center=TRUE,scale=TRUE) */
/* output Z will be standardized */
void standardize(double *Z,int *n,int *k){
  int i,j,ID;
  double *colmean,*colsdev,sumx,sumx2;
  colmean=dvec(*k);
  colsdev=dvec(*k);

  ID = 0;
  for(i=0; i<(*k); i++){
    sumx = 0;
    sumx2 = 0;
    for(j=0; j<(*n); j++){
      colmean[i] = colmean[i] + Z[ID];
      sumx = sumx + Z[ID];
      sumx2 = sumx2 + Z[ID]*Z[ID];
      ID = ID + 1;
    }
    colmean[i] = colmean[i]/(*n);
    colsdev[i] = sqrt((sumx2 - sumx*sumx/(*n))/(*n-1));
  }

  ID = 0;
  for(i=0; i<(*k); i++){
    for(j=0; j<(*n); j++){
      Z[ID] = (Z[ID] - colmean[i])/colsdev[i];
      ID = ID + 1;
    }  
  }
}

/*
algorithm: 
covariance = inverse(precision)
cholesky decomposition of covariance: covariance = LL'
x = mean + Lz
elements in z are iid N(0,1)
*/

/* multivariate normal distribution with mean mu and covariance cov; vec is result. */
void rmvnormal(double *vec, double *mu, double *cov, int *n){
  int i, j, ID, nn,info,incx,incy;
  char *transN="N", *uplo="L";
  double ZERO, ONE,*zvec, *tempcov;
  
  incx = 1;
  incy = 1;
  ZERO = 0;
  ONE = 1;
  info = 0;
  nn = (*n)*(*n);

  zvec=dvec(*n);
  tempcov=dvec(nn);
  GetRNGstate();
  for(i=0; i < *n; i++){
     zvec[i] = rnorm(0, 1);
  } 
  PutRNGstate();
  for(i=0; i < nn; i++){
    tempcov[i] = cov[i];
  }
  /* tempcov = LL'; tempcov contain the L, but the upper triangle is not set to 0*/

  F77_CALL(dpotrf)(uplo,n,tempcov,n,&info FCONE);

  /* printvec(tempcov,nn); */
  /* set the upper triangle to zero; i for column, j for row */
  ID = 0;
  for(i=0; i < *n; i++){
    for(j=0; j < *n; j++){
      if(i > j){
	tempcov[ID] = 0;
      }
      ID = ID + 1;
    }
  }
  /* printvec(tempcov,nn); */
  /* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y */
  /* vec = Lz*/
  F77_CALL(dgemv)(transN,n,n,&ONE,tempcov,n,zvec,&incx, &ZERO,vec, &incy FCONE);  
  dvadd(vec,mu,*n); /* vec = mu +vec */

  Free(zvec);
  Free(tempcov);
} 

/* X is n x p  */
/* Z is n x k  */
/* alpha is p */
/* beta is p x k */
/* sigma2, gamma, accept_gamma are p x 1 */
/* invSigma0 is k+1 x k+1 */
/* invSigmaBeta0 is k+1 */
/* invga0 and invgb0 are scalar */
void bvsNormal(double *X,double *Z,double *alpha,double *beta,double *sigma2,int *gamma,int *accept_gamma,double *prior_ga1,double *invSigma0,double *invSigmaBeta0,double *invga0,double *invgb0,int *n,int *p,int *k){
  char *transN="N", *transT="T";
  double *alphaBeta, *alphaBeta2,*tempab,*beta_m,*beta_V,*beta_invV,*C1Z,*ZtX,*ZtXj,*ZtZ,*invZtZ,*meanX,*meanX_p,*invgamma_b;
  double ONE,ZERO,invgamma_a,sumsq,sumsq_p,prob_gamma,prob_gamma_p,prior_gamma0,prior_gamma1,lhr;
  int *gamma_p,i,j,nk1,k1,k1k1,incx,incy,ID,np;
 
  prior_gamma1 = *prior_ga1;
  prior_gamma0 = 1 - prior_gamma1; 
  k1 = *k + 1;
  nk1 = (*n)*k1;
  np = (*n)*(*p);
  k1k1 = k1*k1;
  
  incx = 1;
  incy = 1;
  ONE = 1;
  ZERO = 0;

  C1Z = dvec(nk1);         /* n x (k+1) cbind(1, Z)*/
  ZtX = dvec((k1)*(*p));  /* (k+1) x p*/
  ZtXj = dvec(k1);
  ZtZ = dvec(k1*k1);       /* (k+1) x (k+1) */
  invZtZ = dvec(k1*k1);     /* (k+1) x (k+1) */
  /* invSigmaBeta0 = dvec(k1); */ /* k+1 */
  meanX = dvec(np);           /* n x p*/
  meanX_p = dvec(np);
  alphaBeta = dvec((*p)*k1);  /* p x (k+1) */
  alphaBeta2 = dvec((*p)*k1); /* p x (k+1) */
  invgamma_b = dvec(*p);         /* p */
  beta_m = dvec(k1);             /* mean beta */
  beta_V = dvec(k1*k1);          /* precision of beta*/
  beta_invV = dvec(k1*k1);       /* variance of beta */
  tempab = dvec(k1);
  gamma_p = ivec(*p);
  /*sigma2 = dvec(*p); */

  /* C1Z = cbind(1, Z) */
  for(i=0; i<nk1; i++){
    if(i < (*n)){
      C1Z[i] = 1;
    }else{
      C1Z[i] = Z[i - (*n)];
    }
  }

  /* compute the posterior mean and variance of beta */
  /*C := alpha*op( A )*op( B ) + beta*C */
  /* ZtZ = t(C1Z) %*% C1Z */
  F77_CALL(dgemm)(transT,transN,&k1,&k1,n,&ONE,C1Z,n,C1Z,n,&ZERO,ZtZ,&k1 FCONE FCONE);
  invsqm2(invZtZ,ZtZ,&k1); /* invZtZ = slove(ZtZ); ZtZ is not changed */
  /*printvec(ZtZ,9);
    printvec(invZtZ,9); */
  /*ZtX = t(C1Z) %*% X; ZtX is k1 by p matrix */
  F77_CALL(dgemm)(transT,transN,&k1,p,n,&ONE,C1Z,n,X,n,&ZERO,ZtX,&k1 FCONE FCONE);
  /* Rprintf("--- ZtX ---\n"); 
     printvec(ZtX,*p*2); */

  /* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y */
  /* invSigmBeta0 = invSigma0 %*% beta0 */
  /* F77_CALL(dgemv)(transN,&k1,&k1,&ONE,invSigma0,&k1,beta0,&incx, &ZERO,invSigmaBeta0, &incy); */ 

  /* cbind(alpha, beta) */
  for(i=0; i<(*p); i++){
    alphaBeta[i] = alpha[i];
  }
  /* scale Beta[p,k] by gamma[p]; each beta[i,] associated gamma[i] */
  ID = 0;
  for(j=0; j<(*k); j++){
    for(i=0; i<(*p); i++){
      /* ID = (*p)*j+i; */
      alphaBeta[(*p)+ID] = beta[ID]*gamma[i];  /* alphaBeta is p by k+1 */
      /* alphaBeta[(*p)+ID] = beta[ID];  because beta is set to 0 if gamma[i]  is 0*/
      ID = ID + 1;
     }
  }

  /* meanx = zz %*% t(alphaBeta)*/
  /* meanX <-  zz[, -1] %*% beta[-1, ] %*% diag.gamma + t(matrix(rep(beta[1, ], n), nrow=p))  */
  /* in R, beta is transformed to (k+1) x p; in C, beta is p x (k+1) */
  F77_CALL(dgemm)(transN,transT,n,p, &k1,&ONE,C1Z,n,alphaBeta,p,&ZERO,meanX,n FCONE FCONE);
  /* Rprintf("--- meanX ---\n");
     printvec(meanX,*p*2); */
  /* posterior parameter for the sigma_j^2 */
  /* meanX = meanX - X; mean X is n x p matrix */
  dvsub(meanX,X,np);
  /* posterior scale parameter b of inverse gamma distribution*/
  ID = 0;
  for(j=0; j<(*p); j++){
    invgamma_b[j] = 0;
    /* t(X_j - Z*beta_j) %*% (X_j - Z*beta_j) */
    for(i=0; i<(*n); i++){
      invgamma_b[j] = invgamma_b[j] + meanX[ID] * meanX[ID];     
      ID = ID + 1;
    }
    invgamma_b[j] = invgamma_b[j]/2 + (*invgb0);
  }

  /* posterior shape parameter a of inverse gamma distribution*/
  invgamma_a = *invga0 + (*n)/2.0;
  
  GetRNGstate();
  /* sampling sigma_j^2 beta_j and gamma_j*/
  ID = 0;
  for(j=0; j<(*p); j++){
    /* Sample sigma2 from its posterior distribution */
    sigma2[j] = rgamma(invgamma_a,1.0/invgamma_b[j]);
    /* posterior variance of beta; gamma*ZtZ*gamma */
    if(gamma[j] == 1){
      dvcopy(beta_V,ZtZ,k1k1);
    }else{
      beta_V[0] = ZtZ[0];
      for(i=1; i<k1k1; i++){
	beta_V[i] = 0;
      }
    }
    
    dvscale(beta_V,k1k1,1.0/sigma2[j]);
    dvadd(beta_V,invSigma0,k1k1);
    invsqm(beta_invV,beta_V,&k1); /* note: beta_V also change */
   
    /* posterior mean of beta */
    ZtXj[0] = ZtX[ID];
    ID = ID + 1;
    for(i=1; i<k1; i++){
      ZtXj[i] = ZtX[ID]*gamma[j];
      ID = ID + 1;
    }
    dvscale(ZtXj,k1,gamma[j]/sigma2[j]);
    dvadd(ZtXj,invSigmaBeta0,k1);
    /* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y */
    /* beta_m = beta_inV %*% ZtXj */
    F77_CALL(dgemv)(transN,&k1,&k1,&ONE,beta_invV,&k1,ZtXj,&incx, &ZERO,beta_m, &incy FCONE);

    /* Generate beta value from multinormal distribution */
    rmvnormal(tempab,beta_m,beta_invV,&k1);
 
    /* Copy tempab to alphaBeta[p,k+1]; */
    for(i=0; i<k1; i++){
      alphaBeta[(*p)*i + j] = tempab[i];
    }

    gamma_p[j] = 1- gamma[j];
  }

  /* copy alphaBeta, which is sampling alpha beta first*/
  for(j=0; j<(*p); j++){
    alpha[j] = alphaBeta[j];
  }

  ID = 0;
  for(j=0; j<(*k); j++){
    for(i=0; i<(*p); i++){
      /* ID = (*p)*j+i; */
      beta[ID] = alphaBeta[(*p)+ID];
      ID = ID + 1;
    }
  }

  /************************** sampling gamma  *************************/
  /* scale Beta[p,k] by gamma[p]; each beta[i,] associated gamma[i] */
  for(j=0; j<(*p); j++){
    alphaBeta2[j] = alphaBeta[j];
  }

  ID = *p;
  for(i=1; i<k1; i++){
    for(j=0; j<(*p); j++){
      /* ID = (*p)*j+i; */
      alphaBeta2[ID] = alphaBeta[ID]*gamma[j]; /* alphaBeta is p by k+1 */
      ID = ID + 1;
    }
  }

  /* meanX <-  zz[, -1] %*% beta[-1, ] %*% diag(cgamma)) + t(matrix(rep(beta[1, ], n), nrow=p))  */
  /* in R, beta is transformed to (k+1) x p; in C, beta is p x (k+1) */
  F77_CALL(dgemm)(transN,transT,n,p, &k1,&ONE,C1Z,n,alphaBeta2,p,&ZERO,meanX,n FCONE FCONE);

  ID = *p;
  for(i=1; i<k1; i++){
    for(j=0; j<(*p); j++){
      /* ID = (*p)*j+i; */
      alphaBeta2[ID] = alphaBeta[ID]*gamma_p[j]; /* alphaBeta is p by k+1 */
      ID = ID + 1;
    }
  }
  
  /* meanX.p <-  zz[, -1] %*% beta[-1, ] %*% diag(c(gamma.p)) + t(matrix(rep(beta[1, ], n), nrow=p))  */
  /* in R, beta is transformed to (k+1) x p; in C, beta is p x (k+1) */
  F77_CALL(dgemm)(transN,transT,n,p, &k1,&ONE,C1Z,n,alphaBeta2,p,&ZERO,meanX_p,n FCONE FCONE);

  /* calculate log acceptance ratio of gamma */
  /* meanX = meanX - X; mean X is n x p matrix */
  dvsub(meanX, X, np);
  dvsub(meanX_p, X, np);
  
  ID = 0;
  for(j=0; j<(*p); j++){
    sumsq = 0;
    sumsq_p = 0;
    for(i=0; i<(*n); i++){
      sumsq = sumsq + meanX[ID] * meanX[ID];
      sumsq_p = sumsq_p + meanX_p[ID] * meanX_p[ID];
      ID = ID + 1;
    }
    sumsq = -0.5*sumsq/sigma2[j];
    sumsq_p = -0.5*sumsq_p/sigma2[j];
    
    prob_gamma = prior_gamma1;
    prob_gamma_p = prior_gamma1;
    if(gamma[j] == 0){
      prob_gamma = prior_gamma0;
    }
    if(gamma_p[j] == 0){
      prob_gamma_p = prior_gamma0;
    }

    lhr = sumsq_p - sumsq + log(prob_gamma_p) - log(prob_gamma);
    if(gamma_p[j] == 0){
      lhr = -1*lhr;
    }
    
    accept_gamma[j] = 0;
    if(runif(0,1) < 1.0/(1.0+exp(-lhr))){
      gamma[j] = 1;
      accept_gamma[j] = 1;
    }else{
      gamma[j] = 0;
    }
  }
  PutRNGstate();
  
  Free(C1Z);         /* n x (k+1) cbind(1, Z)*/
  Free(ZtX);  /* (k+1) x p*/
  Free(ZtXj);
  Free(ZtZ);       /* (k+1) x (k+1) */
  Free(invZtZ);     /* (k+1) x (k+1) */
  Free(meanX);           /* n x p*/
  Free(meanX_p);
  Free(alphaBeta);  /* p x (k+1) */
  Free(alphaBeta2); /* p x (k+1) */
  Free(invgamma_b);         /* p */
  Free(beta_m);             /* mean beta */
  Free(beta_V);          /* variance of beta*/
  Free(beta_invV);       /* precision of beta */
  Free(tempab);
  Free(gamma_p);
}

/* joint sampling gamma and beta */
/* X is n x p  */
/* Z is n x k  */
/* alpha is p */
/* beta is p x k */
/* gamma, accept_gamma are p */
/* invSigma0 is k+1 x k+1 */
/* beta0 is k+1 */
/* lambda_mh is scalar */

void bvsPoisson(int *X,double *Z,double *alpha,double *beta,int *accept_beta,int *gamma,int *accept_gamma,double *prior_ga1,double *invSigma0,double *beta0,double *lambda_mh, int *n,int *p,int *k){
  char *transN="N", *transT="T";
  double *alphaBeta, *alphaBeta2,*alphaBeta_p,*tempab,*tempab_p,*beta_m,*beta_p,*C1Z,*ZtZ,*invZtZ,*ZtBeta,*ZtBeta_p,*Zbeta,*Zbeta_p,*invSigmaBeta0;
  double ONE,ZERO,loglike,loglike_p,prob_gamma,prob_gamma_p,prior_gamma0,prior_gamma1,lhr;
  int *gamma_p,i,j,nk1,k1,incx,incy,ID,np; 

  prior_gamma1 = *prior_ga1;
  prior_gamma0 = 1 - prior_gamma1; 

  k1 = *k + 1;
  nk1 = (*n)*k1;
  np = (*n)*(*p);

  incx = 1;
  incy = 1;
  ONE = 1;
  ZERO = 0;

  C1Z = dvec(nk1);           /* n x (k+1) cbind(1, Z)*/
  ZtZ = dvec(k1*k1);         /* (k+1) x (k+1) */
  invZtZ = dvec(k1*k1);      /* (k+1) x (k+1) */
  invSigmaBeta0 = dvec(k1);  /* k+1 */
  ZtBeta = dvec(np);         /* n x p*/
  ZtBeta_p = dvec(np);
  Zbeta = dvec(*n);         /* n x p*/
  Zbeta_p = dvec(*n); 
  alphaBeta = dvec((*p)*k1);   /* p x (k+1) */
  alphaBeta2 = dvec((*p)*k1);  /* p x (k+1) */
  alphaBeta_p = dvec((*p)*k1); /* p x (k+1) */
  beta_m = dvec(k1);             /* mean beta */
  beta_p = dvec(k1);
  tempab = dvec(k1);
  tempab_p = dvec(k1);
  gamma_p = ivec(*p);

  /* C1Z = cbind(1, Z) */
  for(i=0; i<nk1; i++){
    if(i < (*n)){
      C1Z[i] = 1;
    }else{
      C1Z[i] = Z[i - (*n)];
    }
  }

  /*C := alpha*op( A )*op( B ) + beta*C */
  /* ZtZ = t(C1Z) %*% C1Z */
  F77_CALL(dgemm)(transT,transN,&k1,&k1,n,&ONE,C1Z,n,C1Z,n,&ZERO,ZtZ,&k1 FCONE FCONE);
  invsqm2(invZtZ,ZtZ,&k1); /* invZtZ = slove(ZtZ); ZtZ is not changed */
  dvscale(invZtZ, k1, *lambda_mh); /* var.prop <-  lambda.mh * solve(t(zz) %*% zz); variance for proposal beta */

  /* cbind(alpha, beta) */
  for(i=0; i<(*p); i++){
    alphaBeta[i] = alpha[i];
    alphaBeta2[i] = alpha[i];
  }
  /* scale Beta[p,k] by gamma[p]; each beta[i,] associated gamma[i] */
  ID = 0;
  for(j=0; j<(*k); j++){
    for(i=0; i<(*p); i++){
      /* ID = (*p)*j+i; */
      alphaBeta2[(*p)+ID] = beta[ID]*gamma[i]; /* alphaBeta is p by k+1 */
      alphaBeta[(*p)+ID] = beta[ID]; /* because beta is set to 0 if gamma[i]  is 0*/
      ID = ID + 1;
    }
  }

  /* sampling beta_j and gamma_j */
  /****************** between model move  ***************************/
  GetRNGstate();
  for(j=0; j<(*p); j++){
    /* Generate beta value from multivariate normal distribution */
    /* get beta_j */
    for(i=0; i<k1; i++){
      beta_m[i] = alphaBeta[(*p)*i + j];
    }
    /* generate proposal beta */
    rmvnormal(tempab,beta_m,invZtZ,&k1);
    for(i=0; i<k1; i++){
      alphaBeta_p[(*p)*i + j] = tempab[i];
    }
    gamma_p[j] = 1 - gamma[j];
  }

  /* calculate acceptance ratio for (beta, gamma) jointly */
  /* Ztbeta = zz %*% t(alphaBeta); n x p*/
  /* ZtBeta <-  zz[, -1] %*% beta[-1, ] %*% diag.gamma + t(matrix(rep(beta[1, ], n), nrow=p))  */
  /* in R, beta is transformed to (k+1) x p; in C, beta is p x (k+1) */

  /* calculate ZtBeta using current alphaBeta,which is alphaBeta2 == alpha + (Z)^T beta */
  F77_CALL(dgemm)(transN,transT,n,p, &k1,&ONE,C1Z,n,alphaBeta2,p,&ZERO,ZtBeta,n FCONE FCONE);

  /* calculate ZtBeta_p using proposed alphaBeta,which is alphaBeta2 == alpha + (Z)^T beta */
  for(i=0; i<(*p); i++){
    alphaBeta2[i] = alphaBeta_p[i];
  }
  
  ID = *p;
  for(j=1; j<k1; j++){
    for(i=0; i<(*p); i++){
      alphaBeta2[ID] = alphaBeta_p[ID]*gamma_p[i];
      ID = ID + 1;
    }
  }
  
  F77_CALL(dgemm)(transN,transT,n,p, &k1,&ONE,C1Z,n,alphaBeta2,p,&ZERO,ZtBeta_p,n FCONE FCONE);

  /* calculate joint log-likelihood */
  ID = 0;
  for(j=0; j<(*p); j++){
    loglike = 0;
    loglike_p = 0;
    for(i=0; i<(*n); i++){
      loglike = loglike + X[ID]*ZtBeta[ID] - exp(ZtBeta[ID]);
      loglike_p = loglike_p + X[ID]*ZtBeta_p[ID] - exp(ZtBeta_p[ID]);
      ID = ID + 1;
    }   

    /* Copy tempab to alphaBeta[p,k+1] */
    for(i=0; i<k1; i++){
      tempab[i] = alphaBeta[(*p)*i + j] - beta0[i];
      tempab_p[i] = alphaBeta_p[(*p)*i + j] - beta0[i];
    }
    
   /* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y */
  /* invSigmBeta0 = invSigma0 %*% (beta_j - beta0) */
   F77_CALL(dgemv)(transN,&k1,&k1,&ONE,invSigma0,&k1,tempab,&incx, &ZERO,invSigmaBeta0, &incy FCONE);
   /*  (beta_j - beta0) %*% invSigma0 %*% (beta_j - beta0) */
   loglike = loglike - 0.5 * dvdot(tempab,invSigmaBeta0,k1);

   F77_CALL(dgemv)(transN,&k1,&k1,&ONE,invSigma0,&k1,tempab_p,&incx, &ZERO,invSigmaBeta0, &incy FCONE);  
   loglike_p = loglike_p - 0.5 * dvdot(tempab_p,invSigmaBeta0,k1);

   prob_gamma = prior_gamma1;
   prob_gamma_p = prior_gamma1;
   if(gamma[j] == 0){
     prob_gamma = prior_gamma0;
   }
   if(gamma_p[j] == 0){
     prob_gamma_p = prior_gamma0;
   }

   lhr = loglike_p - loglike + log(prob_gamma_p) - log(prob_gamma); 

   accept_beta[j] = 0;
   accept_gamma[j] = 0;
   if((lhr > 0) | (runif(0,1) < exp(lhr))){ 
     gamma[j] = gamma_p[j];
     alpha[j] = alphaBeta_p[j];
     alphaBeta[j] = alphaBeta_p[j];
     for(i=1; i<k1; i++){
       beta[(*p)*(i-1) + j] = alphaBeta_p[(*p)*i + j];
       alphaBeta[(*p)*i + j] = alphaBeta_p[(*p)*i + j];
     }
     accept_gamma[j] = 1;
     accept_beta[j] = 1;
   }
  }

  /********************* Within model move; do it only when gamma[j] == 1 ***********************/
  for(j=0; j<(*p); j++){
    if(gamma[j] == 1){
      /* Generate beta value from multivariate normal distribution */
      /* get beta_j */
      for(i=0; i<k1; i++){
	beta_m[i] = alphaBeta[(*p)*i + j];
      }
      /* generate proposal beta_p */
      rmvnormal(beta_p,beta_m,invZtZ,&k1);

      /* calculate acceptance ratio for beta */
      /* Ztbeta = zz %*% alpha_beta; n x 1*/
      /* calculate Ztbeta using current alpha_beta,which is alphaBeta2 == alpha + (Z)^T beta */
      F77_CALL(dgemv)(transN,n,&k1,&ONE,C1Z,n,beta_m,&incx,&ZERO,Zbeta,&incy FCONE);
      
      /* calculate ZtBeta_p using proposed alphaBeta,which is alphaBeta2 == alpha + (Z)^T beta */
      F77_CALL(dgemv)(transN,n,&k1,&ONE,C1Z,n,beta_p,&incx,&ZERO,Zbeta_p,&incy FCONE);

      /* calculate joint log-likelihood;  note X is n x p */
      loglike = 0;
      loglike_p = 0;
      for(i=0; i<(*n); i++){
	ID = (*n)*j + i;
	loglike = loglike + X[ID]*Zbeta[i] - exp(Zbeta[i]);
	loglike_p = loglike_p + X[ID]*Zbeta_p[i] - exp(Zbeta_p[i]);
      }   

      /* Copy tempab to alphaBeta[p,k+1] */
      for(i=0; i<k1; i++){
	tempab[i] = beta_m[i] - beta0[i];
	tempab_p[i] = beta_p[i] - beta0[i];
      }
    
      /* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y */
      /* invSigmBeta0 = invSigma0 %*% (beta_j - beta0) */
      F77_CALL(dgemv)(transN,&k1,&k1,&ONE,invSigma0,&k1,tempab,&incx, &ZERO,invSigmaBeta0, &incy FCONE);
      /*  (beta_j - beta0) %*% invSigma0 %*% (beta_j - beta0) */
      loglike = loglike - 0.5 * dvdot(tempab,invSigmaBeta0,k1);

      F77_CALL(dgemv)(transN,&k1,&k1,&ONE,invSigma0,&k1,tempab_p,&incx, &ZERO,invSigmaBeta0, &incy FCONE);  
      loglike_p = loglike_p - 0.5 * dvdot(tempab_p,invSigmaBeta0,k1);
      /* don't need to calculate log probability of gamma[j] because the current and proposed gamma[j] are 1 */
      lhr = loglike_p - loglike; 

      accept_beta[j] = 0;
      if((lhr > 0) | (runif(0,1) < exp(lhr))){ 
	alpha[j] = beta_p[0];
	for(i=1; i<k1; i++){
	  beta[(*p)*(i-1) + j] = beta_p[i];
	}
	accept_beta[j] = 1;
      }
    }
  }
 
  PutRNGstate();
    
  /* ********************************************************************/
  
  Free(C1Z);                   /* n x (k+1) cbind(1, Z)*/
  Free(ZtZ);                  /* (k+1) x (k+1) */
  Free(invZtZ);              /* (k+1) x (k+1) */
  Free(invSigmaBeta0);      /* k+1 */
  Free(ZtBeta);            /* n x p*/
  Free(ZtBeta_p);         /* n x p*/
  Free(Zbeta);            /* p*/
  Free(Zbeta_p);         /* p*/
  Free(alphaBeta);       /* p x (k+1) */
  Free(alphaBeta2);     /* p x (k+1) */
  Free(alphaBeta_p);   /* p x (k+1) */
  Free(beta_m);       /* mean beta */
  Free(beta_p);
  Free(tempab);
  Free(tempab_p);
  Free(gamma_p);
}


/* jointly sample alphaBeta and gamma */
/* X is n x p  */
/* Z is n x k  */
/* alpha is p */
/* beta is p x k */
/* gamma, accept_gamma are p */
/* invSigma0 is k+1 x k+1 */
/* beta0 is k+1 */
/* lambda_mh is scalar */
/* accept_beta is p */

void bvsBinom(int *X,double *Z,double *alpha,double *beta,int *accept_beta,int *gamma,int *accept_gamma,double *prior_ga1,double *invSigma0,double *beta0,double *lambda_mh, int *n,int *p,int *k){
  char *transN="N", *transT="T";
  double *alphaBeta, *alphaBeta2,*alphaBeta_p,*tempab,*tempab_p,*beta_m,*beta_p,*C1Z,*ZtZ,*invZtZ,*ZtBeta,*ZtBeta_p,*Zbeta,*Zbeta_p,*invSigmaBeta0;
  double ONE,ZERO,loglike,loglike_p,prob_gamma,prob_gamma_p,prior_gamma0,prior_gamma1,lhr;
  int *gamma_p,i,j,nk1,k1,incx,incy,ID,np;

  prior_gamma1 = *prior_ga1;
  prior_gamma0 = 1 - prior_gamma1; 

  k1 = *k + 1;
  nk1 = (*n)*k1;
  np = (*n)*(*p);

  incx = 1;
  incy = 1;
  ONE = 1;
  ZERO = 0;

  C1Z = dvec(nk1);           /* n x (k+1) cbind(1, Z)*/
  ZtZ = dvec(k1*k1);         /* (k+1) x (k+1) */
  invZtZ = dvec(k1*k1);      /* (k+1) x (k+1) */
  invSigmaBeta0 = dvec(k1);  /* k+1 */
  ZtBeta = dvec(np);         /* n x p*/
  ZtBeta_p = dvec(np);
  Zbeta = dvec(*n);         /* n x p*/
  Zbeta_p = dvec(*n); 
  alphaBeta = dvec((*p)*k1);   /* p x (k+1) */
  alphaBeta2 = dvec((*p)*k1);  /* p x (k+1) */
  alphaBeta_p = dvec((*p)*k1); /* p x (k+1) */
  beta_m = dvec(k1);             /* mean beta */
  beta_p = dvec(k1);
  tempab = dvec(k1);
  tempab_p = dvec(k1);
  gamma_p = ivec(*p);

  /* C1Z = cbind(1, Z) */
  for(i=0; i<nk1; i++){
    if(i < (*n)){
      C1Z[i] = 1;
    }else{
      C1Z[i] = Z[i - (*n)];
    }
  }

  /*C := alpha*op( A )*op( B ) + beta*C */
  /* ZtZ = t(C1Z) %*% C1Z */
  F77_CALL(dgemm)(transT,transN,&k1,&k1,n,&ONE,C1Z,n,C1Z,n,&ZERO,ZtZ,&k1 FCONE FCONE);
  invsqm2(invZtZ,ZtZ,&k1); /* invZtZ = slove(ZtZ); ZtZ is not changed */
  dvscale(invZtZ, k1, *lambda_mh); /* var.prop <-  lambda.mh * solve(t(zz) %*% zz); variance for proposal beta */

    /* cbind(alpha, beta) */
  for(i=0; i<(*p); i++){
    alphaBeta[i] = alpha[i];
    alphaBeta2[i] = alpha[i];
  }
  /* scale Beta[p,k] by gamma[p]; each beta[i,] associated gamma[i] */
  ID = 0;
  for(j=0; j<(*k); j++){
    for(i=0; i<(*p); i++){
      /* ID = (*p)*j+i; */
      alphaBeta2[(*p)+ID] = beta[ID]*gamma[i];
      alphaBeta[(*p)+ID] = beta[ID]; 
      ID = ID + 1;
    }
  }

  /* sampling beta_j and gamma_j */
  /****************** between model move  ***************************/
  GetRNGstate();
  for(j=0; j<(*p); j++){
    /* Generate beta value from multivariate normal distribution */
    /* get beta_j */
    for(i=0; i<k1; i++){
      beta_m[i] = alphaBeta[(*p)*i + j];
    }
    /* generate proposal beta; alphaBeta is p x (k+1) */
    rmvnormal(tempab,beta_m,invZtZ,&k1);
    for(i=0; i<k1; i++){
      alphaBeta_p[(*p)*i + j] = tempab[i];
    }
    gamma_p[j] = 1 - gamma[j];
  }

  /* calculate acceptance ratio for (beta, gamma) jointly */
  /* Ztbeta = zz %*% t(alphaBeta); n x p*/
  /* ZtBeta <-  zz[, -1] %*% beta[-1, ] %*% diag.gamma + t(matrix(rep(beta[1, ], n), nrow=p))  */
  /* in R, beta is transformed to (k+1) x p; in C, beta is p x (k+1) */

  /* calculate ZtBeta using current alphaBeta,which is alphaBeta2 == alpha + (Z)^T beta */
  F77_CALL(dgemm)(transN,transT,n,p, &k1,&ONE,C1Z,n,alphaBeta2,p,&ZERO,ZtBeta,n FCONE FCONE);

  /* calculate ZtBeta_p using proposed alphaBeta,which is alphaBeta2 == alpha + (Z)^T beta */
  for(i=0; i<(*p); i++){
    alphaBeta2[i] = alphaBeta_p[i];
  }
  
  ID = *p;
  for(j=1; j<k1; j++){
    for(i=0; i<(*p); i++){
      alphaBeta2[ID] = alphaBeta_p[ID]*gamma_p[i];
      ID = ID + 1;
    }
  }
  
  F77_CALL(dgemm)(transN,transT,n,p, &k1,&ONE,C1Z,n,alphaBeta2,p,&ZERO,ZtBeta_p,n FCONE FCONE);

  /*********** calculate joint log likelihood *************/
  ID = 0;
  for(j=0; j<(*p); j++){
    loglike = 0;
    loglike_p = 0;
    for(i=0; i<(*n); i++){
      if(X[ID]==1){
	loglike = loglike + logpnull(ZtBeta[ID]);
	loglike_p = loglike_p + logpnull(ZtBeta_p[ID]);
      }else{
	loglike = loglike + logqnull(ZtBeta[ID]);
	loglike_p = loglike_p + logqnull(ZtBeta_p[ID]);
      }
      ID = ID + 1;
    }   

    /* Copy tempab to alphaBeta[p,k+1] */
    for(i=0; i<k1; i++){
      tempab[i] = alphaBeta[(*p)*i + j] - beta0[i];
      tempab_p[i] = alphaBeta_p[(*p)*i + j] - beta0[i];
    }
    
   /* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y */
  /* invSigmBeta0 = invSigma0 %*% (beta_j - beta0) */
   F77_CALL(dgemv)(transN,&k1,&k1,&ONE,invSigma0,&k1,tempab,&incx, &ZERO,invSigmaBeta0, &incy FCONE);  
   loglike = loglike - 0.5 * dvdot(tempab,invSigmaBeta0,k1);

   F77_CALL(dgemv)(transN,&k1,&k1,&ONE,invSigma0,&k1,tempab_p,&incx, &ZERO,invSigmaBeta0, &incy FCONE);  
   loglike_p = loglike_p - 0.5 * dvdot(tempab_p,invSigmaBeta0,k1);

   prob_gamma = prior_gamma1;
   prob_gamma_p = prior_gamma1;
   if(gamma[j] == 0){
     prob_gamma = prior_gamma0;
   }
   if(gamma_p[j] == 0){
     prob_gamma_p = prior_gamma0;
   }

   lhr = loglike_p - loglike + log(prob_gamma_p) - log(prob_gamma); 

   accept_beta[j] = 0;
   accept_gamma[j] = 0;
   if((lhr > 0) | (runif(0,1) < exp(lhr))){ 
     gamma[j] = gamma_p[j];
     alpha[j] = alphaBeta_p[j];
     alphaBeta[j] = alphaBeta_p[j];
     for(i=1; i<k1; i++){
       beta[(*p)*(i-1) + j] = alphaBeta_p[(*p)*i + j];
       alphaBeta[(*p)*i + j] = alphaBeta_p[(*p)*i + j];
     }
     accept_gamma[j] = 1;
     accept_beta[j] = 1;
   }
  }

  /** Within model move; do it only when current gamma[j] == 1, the proposed new gamma[j] is also 1 **/
  for(j=0; j<(*p); j++){
    if(gamma[j] == 1){
      /* Generate beta value from multivariate normal distribution */
      /* get beta_j */
      for(i=0; i<k1; i++){
	beta_m[i] = alphaBeta[(*p)*i + j];
      }
      /* generate proposal beta */
      rmvnormal(beta_p,beta_m,invZtZ,&k1);

      /* calculate acceptance ratio for beta */
      /* Ztbeta = zz %*% alpha_beta; n x 1*/
      /* calculate Ztbeta using current alpha_beta,which is alphaBeta2 == alpha + (Z)^T beta */
      F77_CALL(dgemv)(transN,n,&k1,&ONE,C1Z,n,beta_m,&incx,&ZERO,Zbeta,&incy FCONE);
      
      /* calculate ZtBeta_p using proposed alphaBeta,which is alphaBeta2 == alpha + (Z)^T beta */
      F77_CALL(dgemv)(transN,n,&k1,&ONE,C1Z,n,beta_p,&incx,&ZERO,Zbeta_p,&incy FCONE);

      /* calculate joint log-likelihood */
      loglike = 0;
      loglike_p = 0;
      for(i=0; i<(*n); i++){
	ID = (*n)*j + i;
	if(X[ID]==1){
	  loglike = loglike + logpnull(Zbeta[i]);
	  loglike_p = loglike_p + logpnull(Zbeta_p[i]);
	}else{
	  loglike = loglike + logqnull(Zbeta[i]);
	  loglike_p = loglike_p + logqnull(Zbeta_p[i]);
	}
      }   

      /* Copy tempab to alphaBeta[p,k+1] */
      for(i=0; i<k1; i++){
	tempab[i] = beta_m[i] - beta0[i];
	tempab_p[i] = beta_p[i] - beta0[i];
      }
    
      /* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y */
      /* invSigmBeta0 = invSigma0 %*% (beta_j - beta0) */
      F77_CALL(dgemv)(transN,&k1,&k1,&ONE,invSigma0,&k1,tempab,&incx, &ZERO,invSigmaBeta0, &incy FCONE);
      /*  (beta_j - beta0) %*% invSigma0 %*% (beta_j - beta0) */
      loglike = loglike - 0.5 * dvdot(tempab,invSigmaBeta0,k1);

      F77_CALL(dgemv)(transN,&k1,&k1,&ONE,invSigma0,&k1,tempab_p,&incx, &ZERO,invSigmaBeta0, &incy FCONE);  
      loglike_p = loglike_p - 0.5 * dvdot(tempab_p,invSigmaBeta0,k1);

      /* Don't need to calculate the log probability of gamma[j] because both current and proposed gamma[j] are 1. */
      lhr = loglike_p - loglike;

      accept_beta[j] = 0; 
      if((lhr > 0) | (runif(0,1) < exp(lhr))){ 
	alpha[j] = beta_p[0];
	for(i=1; i<k1; i++){
	  beta[(*p)*(i-1) + j] = beta_p[i];
	}
	accept_beta[j] = 1;
      }
    }
  }
 
  PutRNGstate();

  Free(C1Z);                   /* n x (k+1) cbind(1, Z)*/
  Free(ZtZ);                  /* (k+1) x (k+1) */
  Free(invZtZ);              /* (k+1) x (k+1) */
  Free(invSigmaBeta0);      /* k+1 */
  Free(ZtBeta);            /* n x p*/
  Free(ZtBeta_p);         /* n x p*/
  Free(Zbeta);            /* p*/
  Free(Zbeta_p);         /* p*/
  Free(alphaBeta);       /* p x (k+1) */
  Free(alphaBeta2);     /* p x (k+1) */
  Free(alphaBeta_p);   /* p x (k+1) */
  Free(beta_m);       /* mean beta */
  Free(beta_p);
  Free(tempab);
  Free(tempab_p);
  Free(gamma_p);
}

/*if tryz is accepted, accept = accept + 1; function for iClusterBayes */
void metroMix6d(double *zi,dataType *dt1,dataType *dt2,dataType *dt3,dataType *dt4,dataType *dt5,dataType *dt6,double *sigma,
	      int *ndt,int *accept){
  double dif, *tryz;
  int i,k,incx, incy;

  incx = 1;
  incy = 1;

  k = *(dt1->k);
  tryz = dvec(k);

  for(i=0; i<k;i++){
    tryz[i] = zi[i] + rnorm(0, *sigma);
  }
  /* case 4 is not used in iClusterBayes */
  dif = 0;
  if((*ndt) > 0){
    switch(*(dt1->type)){
    case 1 : dif += logNorm(tryz,dt1->alpha,dt1->beta,dt1->sigma2,dt1->con,dt1->p,dt1->k) 
	- logNorm(zi,dt1->alpha,dt1->beta,dt1->sigma2,dt1->con,dt1->p,dt1->k); break;
      
    case 2 : dif += logBinom(tryz,dt1->alpha,dt1->beta,dt1->cat,dt1->p,dt1->k) 
	- logBinom(zi,dt1->alpha,dt1->beta,dt1->cat,dt1->p,dt1->k); break;

    case 3 : dif += logPoisson(tryz,dt1->alpha,dt1->beta,dt1->cat,dt1->p,dt1->k) 
	- logPoisson(zi,dt1->alpha,dt1->beta,dt1->cat,dt1->p,dt1->k); break;
    }
  }

  if((*ndt) > 1){
    switch(*(dt2->type)){
    case 1 : dif += logNorm(tryz,dt2->alpha,dt2->beta,dt2->sigma2,dt2->con,dt2->p,dt2->k) 
	- logNorm(zi,dt2->alpha,dt2->beta,dt2->sigma2,dt2->con,dt2->p,dt2->k); break;
      
    case 2 : dif += logBinom(tryz,dt2->alpha,dt2->beta,dt2->cat,dt2->p,dt2->k) 
	- logBinom(zi,dt2->alpha,dt2->beta,dt2->cat,dt2->p,dt2->k); break;

    case 3 : dif += logPoisson(tryz,dt2->alpha,dt2->beta,dt2->cat,dt2->p,dt2->k) 
	- logPoisson(zi,dt2->alpha,dt2->beta,dt2->cat,dt2->p,dt2->k); break;
    }
   }
 
  if((*ndt) > 2){
    switch(*(dt3->type)){
    case 1 : dif += logNorm(tryz,dt3->alpha,dt3->beta,dt3->sigma2,dt3->con,dt3->p,dt3->k) 
	- logNorm(zi,dt3->alpha,dt3->beta,dt3->sigma2,dt3->con,dt3->p,dt3->k); break;
      
    case 2 : dif += logBinom(tryz,dt3->alpha,dt3->beta,dt3->cat,dt3->p,dt3->k) 
	- logBinom(zi,dt3->alpha,dt3->beta,dt3->cat,dt3->p,dt3->k); break;
      
    case 3 : dif += logPoisson(tryz,dt3->alpha,dt3->beta,dt3->cat,dt3->p,dt3->k) 
	- logPoisson(zi,dt3->alpha,dt3->beta,dt3->cat,dt3->p,dt3->k); break;
    }
   }
  
  if((*ndt) > 3){
    switch(*(dt4->type)){
    case 1 : dif += logNorm(tryz,dt4->alpha,dt4->beta,dt4->sigma2,dt4->con,dt4->p,dt4->k) 
	- logNorm(zi,dt4->alpha,dt4->beta,dt4->sigma2,dt4->con,dt4->p,dt4->k); break;
      
    case 2 : dif += logBinom(tryz,dt4->alpha,dt4->beta,dt4->cat,dt4->p,dt4->k) 
	- logBinom(zi,dt4->alpha,dt4->beta,dt4->cat,dt4->p,dt4->k); break;
      
    case 3 : dif += logPoisson(tryz,dt4->alpha,dt4->beta,dt4->cat,dt4->p,dt4->k) 
	- logPoisson(zi,dt4->alpha,dt4->beta,dt4->cat,dt4->p,dt4->k); break;
    }
  }

  if((*ndt) > 4){
    switch(*(dt5->type)){
    case 1 : dif += logNorm(tryz,dt5->alpha,dt5->beta,dt5->sigma2,dt5->con,dt5->p,dt5->k) 
	- logNorm(zi,dt5->alpha,dt5->beta,dt5->sigma2,dt5->con,dt5->p,dt5->k); break;
      
    case 2 : dif += logBinom(tryz,dt5->alpha,dt5->beta,dt5->cat,dt5->p,dt5->k) 
	- logBinom(zi,dt5->alpha,dt5->beta,dt5->cat,dt5->p,dt5->k); break;
      
    case 3 : dif += logPoisson(tryz,dt5->alpha,dt5->beta,dt5->cat,dt5->p,dt5->k) 
	- logPoisson(zi,dt5->alpha,dt5->beta,dt5->cat,dt5->p,dt5->k); break;
    }
  }

  if((*ndt) > 5){
    switch(*(dt6->type)){
    case 1 : dif += logNorm(tryz,dt6->alpha,dt6->beta,dt6->sigma2,dt6->con,dt6->p,dt6->k) 
	- logNorm(zi,dt6->alpha,dt6->beta,dt6->sigma2,dt6->con,dt6->p,dt6->k); break;
      
    case 2 : dif += logBinom(tryz,dt6->alpha,dt6->beta,dt6->cat,dt6->p,dt6->k) 
	- logBinom(zi,dt6->alpha,dt6->beta,dt6->cat,dt6->p,dt6->k); break;
      
    case 3 : dif += logPoisson(tryz,dt6->alpha,dt6->beta,dt6->cat,dt6->p,dt6->k) 
	- logPoisson(zi,dt6->alpha,dt6->beta,dt6->cat,dt6->p,dt6->k); break;
    }
  }
    
  dif = dif - 0.5*F77_CALL(ddot)(&k,tryz,&incx,tryz,&incy) + 0.5*F77_CALL(ddot)(&k,zi,&incx,zi,&incy);
  
  if((dif > 0) | (runif(0,1) < exp(dif))){
    dvcopy(zi, tryz, k);
    *accept = *accept + 1;
  }

  Free(tryz); 
}

/* This function is modified from mcmcMix in the iClusterPlus.c; this fun can handle 6 data sets */
/* lastZ is changed after running the function */
/* meanZ is the mean of latent variable Z, lastZ is the last draw of Z; when ndraw=1,meanZ==lastZ,
   nkbd[0] is n; nkbd[1] is k; nkbd[2] is burnin; nkbd[3] is draw; nkbd[4] is ndt;
   n: number of subjects; k: number of clusters; sdev - used to control random walk variation; 
   ndt: number of data types; typ; data type, 1 is normal, 2 is binary, 3 poisson,(4 multinomial--not modeled in iClusterBayes); 
   a0-5 is p dimensional,b0-5[p][k],con0-5[n][p] or cat0-5[n][p], sigma0-5[p] 
   con0-3: continuous response; cat0-3: categorical response or Poisson count; 
   accept: n x 1 vector (should be initialized as 0) to track if the proposal latent variable is accepted 
*/
void mcmcMix6d(double *meanZ,double *lastZ,int *nkbd,double *sdev, 
	       int *ty0,int *p0,double *a0,double *b0,double *con0,int *cat0,double *sigma0,
	       int *ty1,int *p1,double *a1,double *b1,double *con1,int *cat1,double *sigma1,
	       int *ty2,int *p2,double *a2,double *b2,double *con2,int *cat2,double *sigma2,
	       int *ty3,int *p3,double *a3,double *b3,double *con3,int *cat3,double *sigma3,
	       int *ty4,int *p4,double *a4,double *b4,double *con4,int *cat4,double *sigma4,
	       int *ty5,int *p5,double *a5,double *b5,double *con5,int *cat5,double *sigma5,
	       int *accept){

  double **newz, **tempz;
  double **conList[6];
  int **catList[6];
  int n,k,burnin,draw,ndt,i,j,t,ID;
  int conp[6] = {1,1,1,1,1,1}; 
  int conn[6] = {1,1,1,1,1,1};
  int catp[6] = {1,1,1,1,1,1};
  int catn[6] = {1,1,1,1,1,1};
  dataType *dt[6];

  n = nkbd[0];
  k = nkbd[1];
  burnin = nkbd[2];
  draw = nkbd[3];
  ndt = nkbd[4];
  
  /* decide the size of allocated matrix */
  if((*ty0) == 1){
    conp[0] = *p0;
    conn[0] = n;
  }else{
    catp[0] = *p0;
    catn[0] = n;
  }

  if((*ty1) == 1){
    conp[1] = *p1;
    conn[1] = n;
  }else{
    catp[1] = *p1;
    catn[1] = n;
  }

  if((*ty2) == 1){
    conp[2] = *p2;
    conn[2] = n;
  }else{
    catp[2] = *p2;
    catn[2] = n;
  }

  if((*ty3) == 1){
    conp[3] = *p3;
    conn[3] = n;
  }else{
    catp[3] = *p3;
    catn[3] = n;
  }

  if((*ty4) == 1){
    conp[4] = *p4;
    conn[4] = n;
  }else{
    catp[4] = *p4;
    catn[4] = n;
  }

  if((*ty5) == 1){
    conp[5] = *p5;
    conn[5] = n;
  }else{
    catp[5] = *p5;
    catn[5] = n;
  }
  
  /*  nrow = (n)*(draw); */
  newz = drowm(n, k);
  tempz = drowm(n, k);

  for(i=0; i<6; i++){
    dt[i] = (dataType *)malloc(sizeof(dataType));
    if(dt[i] == NULL){
      error("Error: cannot allocate memory for dt[]\n");
      /*    exit(1); */
    }
    conList[i] = drowm(conn[i],conp[i]);
    catList[i] = irowm(catn[i],catp[i]);
  }
  /* assign the values column by column */
  ID = 0;
  for(j=0; j<(k); j++){
    for(i=0; i<(n); i++){
      tempz[i][j] = lastZ[ID]; 
      ID = ID + 1;
    }
  }

  if((*ty0) == 1){  /* 1 is normal data */
    dvtom(conList[0],con0,n,*p0); /* transform response data to n by p0 matrix */
  }else{             /* 2,3,4 is integer data */
    ivtom(catList[0],cat0,n,*p0);  /* transform response data to n by p0 matrix */
  }   /*conList[0][0] and catList[0][0] is px1 vector; dt[0]->con or dt[0]->cat will be replaced in MCMC */
  fillData3(dt[0],a0,b0,conList[0][0],catList[0][0],sigma0,p0,&k,ty0); 
  
  /* for the second data set */
  if((ndt) > 1){
    if((*ty1) == 1){
      dvtom(conList[1],con1,n,*p1);
    }else{
      ivtom(catList[1],cat1,n,*p1);
    }
    fillData3(dt[1],a1,b1,conList[1][0],catList[1][0],sigma1,p1,&k,ty1);
  }
  /* for the third data set */  
  if((ndt) > 2){
    if((*ty2) == 1){
      dvtom(conList[2],con2,n,*p2);
    }else{
      ivtom(catList[2],cat2,n,*p2);
    }
    fillData3(dt[2],a2,b2,conList[2][0],catList[2][0],sigma2,p2,&k,ty2);
  }
  /* for the four data set */
  if((ndt) > 3){
    if((*ty3) == 1){
      dvtom(conList[3],con3,n,*p3);
    }else{
      ivtom(catList[3],cat3,n,*p3);
    }
    fillData3(dt[3],a3,b3,conList[3][0],catList[3][0],sigma3,p3,&k,ty3);
  }
  /* for the five data set */
  if((ndt) > 4){
    if((*ty4) == 1){
      dvtom(conList[4],con4,n,*p4);
    }else{
      ivtom(catList[4],cat4,n,*p4);
    }
    fillData3(dt[4],a4,b4,conList[4][0],catList[4][0],sigma4,p4,&k,ty4);
  }
  /* for the six data set */
  if((ndt) > 5){
    if((*ty5) == 1){
      dvtom(conList[5],con5,n,*p5);
    }else{
      ivtom(catList[5],cat5,n,*p5);
    }
    fillData3(dt[5],a5,b5,conList[5][0],catList[5][0],sigma5,p5,&k,ty5);
  }

  /* Important to use GetRNGstate */
  GetRNGstate();
  /* MCMC sampling */
  /* burn.in + the first draw*/
  for(i=0; i<(burnin); i++){
    for(j=0; j<(n); j++){
      /* prepare the data from 1 to ndt */
      for(t=0; t<ndt; t++){
	if(*(dt[t]->type) == 1){
	  dt[t]->con = conList[t][j];
	}else{
	  dt[t]->cat = catList[t][j];
	}
      }
      metroMix6d(tempz[j],dt[0],dt[1],dt[2],dt[3],dt[4],dt[5],sdev,&ndt,&accept[j]);
    }
  }
  
  /* sampling */
  for(i=0; i<(draw); i++){
    for(j=0; j<(n); j++){
      /* prepare the data from 1 to ndt */
      for(t=0; t<ndt; t++){
	if(*(dt[t]->type) == 1){
	  dt[t]->con = conList[t][j];
	}else{
	  dt[t]->cat = catList[t][j];
	}
      }
      metroMix6d(tempz[j],dt[0],dt[1],dt[2],dt[3],dt[4],dt[5],sdev,&ndt,&accept[j]); 
    }
    dmadd(newz,tempz,n,k);
  }
 
  PutRNGstate();
  if((draw)>1){
    dmscale(newz,n,k,1.0/(draw)); /* average of Z */
  }
  dmtov(meanZ,newz,n,k);
  dmtov(lastZ,tempz,n,k);
  dmfree(newz, n);
  dmfree(tempz, n);
  
  for(i=0; i<6; i++){
    dmfree(conList[i],conn[i]);
    imfree(catList[i],catn[i]);
    free(dt[i]);      
  }

}

/* if the data type is multinomial, beta = p  x (k x C) is the format in R; 
   beta must be transposed using t(beta) before passing to this function */
/* lastZ is changed after running the function */
/* meanZ is the mean of latent variable Z, lastZ is the last draw of Z;
   n: number of subjects; k: number of clusters; sdev - used to control random walk variation; 
   ndt: number of data types; typ; data type, 1 is normal, 2 is binary, 3 poisson, 4 multinomial; 
   a0-5 is p dimensional,b0-5[p][k],con0-5[n][p] or cat0-5[n][p], sigma0-5[p] 
   con0-3: continuous response; cat0-3: categorical response; 
   acceptZ: n x 1 vector (should be initialized as 0) to track if the proposal latent variable is accepted 
   nkZbd =c(n,k,Zburnin,Zdraw,ndt,thin)
   thin is used to control the dense of sumMeanZ and suma0-3,sumb0-3; if betaDraw % thin == 0, add meanZ,a0-3,b0-3,respectively.
   pg is ndt-dimensional vector for priors of gamma
   invSigma0 is k+1 x k+1 
   invSigmaBeta0 is k+1 
   invgab0_lbd for inverse gamma prior a0, b0, and the scaler for the covariance matrix of beta proposal 
   gamma0-5,acceptBeta0-5,acceptGamma0-5 are all p-dimensional vector,
*/
void mcmcBayes(double *meanZ,double *lastZ,int *acceptZ,int *nkZbd,double *sdev,
	       int *typc0,double *a0,double *b0,double *con0,int *cat0,double *sigma0,
	       int *typc1,double *a1,double *b1,double *con1,int *cat1,double *sigma1,
	       int *typc2,double *a2,double *b2,double *con2,int *cat2,double *sigma2,
	       int *typc3,double *a3,double *b3,double *con3,int *cat3,double *sigma3,
	       int *typc4,double *a4,double *b4,double *con4,int *cat4,double *sigma4,
	       int *typc5,double *a5,double *b5,double *con5,int *cat5,double *sigma5,
	       int *betaBurninDraw, double *beta0,double *invSigma0,double *invSigmaBeta0,double *invgab0_lbd, double*pg,
	       int *gamma0,int *acceptBeta0,int *acceptGamma0,int *gamma1,int *acceptBeta1,int *acceptGamma1,
	       int *gamma2,int *acceptBeta2,int *acceptGamma2,int *gamma3,int *acceptBeta3,int *acceptGamma3,
	       int *gamma4,int *acceptBeta4,int *acceptGamma4,int *gamma5,int *acceptBeta5,int *acceptGamma5){

  int i,j,h,n,k,ID,Zburnin,Zdraw,ndt,nk,thin,nthin,remainder,ty0,p0,ty1,p1,ty2,p2,ty3,p3,ty4,p4,ty5,p5,betaBurnin,betaDraw,p0k,p1k,p2k,p3k,p4k,p5k;
  double invga0,invgb0,mhLambda,scale;
  double *sumMeanZ,*suma0,*sumb0,*suma1,*sumb1,*suma2,*sumb2,*suma3,*sumb3,*suma4,*sumb4,*suma5,*sumb5,*sumsig0,*sumsig1,*sumsig2,*sumsig3,*sumsig4,*sumsig5,*gb0,*gb1,*gb2,*gb3,*gb4,*gb5;
  int *sumGamma0,*sumAcceptGamma0,*sumAcceptBeta0,*sumGamma1,*sumAcceptGamma1,*sumAcceptBeta1,*sumGamma2,*sumAcceptGamma2,*sumAcceptBeta2;
  int *sumGamma3,*sumAcceptGamma3,*sumAcceptBeta3,*sumGamma4,*sumAcceptGamma4,*sumAcceptBeta4,*sumGamma5,*sumAcceptGamma5,*sumAcceptBeta5;
  /* int class0=1,nclass0=1,class1=1,nclass1=1,class2=1,nclass2=1,class3=1,nclass3=1,class4=1,nclass4=1,class5=1,nclass5=1; */
  /* class0-5 and nclass0-5 are used for multi-categorical data, which are not moldeled by this bayesian method */
  /* The reason of keeping class0-5 and nclass0-5 is just in order to use previous iClusterPlus code */
  
  n = nkZbd[0];
  k = nkZbd[1];
  Zburnin = nkZbd[2];
  Zdraw = nkZbd[3];
  ndt = nkZbd[4];
  thin = nkZbd[5];
  nk = n*k;
  nthin = 0;

  ty0 = typc0[0];
  p0 = typc0[1];
  
  ty1 = typc1[0];
  p1 = typc1[1];

  ty2 = typc2[0];
  p2 = typc2[1];

  ty3 = typc3[0];
  p3 = typc3[1];

  ty4 = typc4[0];
  p4 = typc4[1];
  
  ty5 = typc5[0];
  p5 = typc5[1];

  /* if ndt < 3; the following operation such as dvadd(sumb3,b3,p3k) will casue segmentation fault because the size of b3 < p3k;*/
  /* set p*k to 1 to avoid segmentation fault */
  p0k = 1;
  p1k = 1;
  p2k = 1;
  p3k = 1;
  p4k = 1;
  p5k = 1;

  if(ndt > 0){
    p0k = p0*k;
    gb0 = dvec(p0k);
    dvcopy(gb0,b0,p0k);
  }
  if(ndt > 1){
    p1k = p1*k;
    gb1 = dvec(p1k);
    dvcopy(gb1,b1,p1k);
  }
  if(ndt > 2){
    p2k = p2*k;
    gb2 = dvec(p2k);
    dvcopy(gb2,b2,p2k);
  }
  if(ndt > 3){
    p3k = p3*k;
    gb3 = dvec(p3k);
    dvcopy(gb3,b3,p3k);
  }
  if(ndt > 4){
    p4k = p4*k;
    gb4 = dvec(p4k);
    dvcopy(gb4,b4,p4k);
  }
  if(ndt > 5){
    p5k = p5*k;
    gb5 = dvec(p5k);
    dvcopy(gb5,b5,p5k);
  }

  betaBurnin = betaBurninDraw[0];
  betaDraw = betaBurninDraw[1];

  invga0 = invgab0_lbd[0];
  invgb0 = invgab0_lbd[1];
  mhLambda = invgab0_lbd[2];

  /* if data type is not normal, length of p and sigma are not the same, which call segmentation fault */
  /* 1 is normal case; 2 is binomial; 3 is poisson; length(sigma2) is 1 */
  if(ty0 == 1){
    sumsig0 = dvec(p0);
  }else{
    sumsig0 = dvec(1);
  }
  suma0 = dvec(p0);
  sumb0 = dvec(p0k);
  sumGamma0 = ivec(p0);
  sumAcceptGamma0 = ivec(p0);
  sumAcceptBeta0 = ivec(p0);

  if(ty1 == 1){
    sumsig1 = dvec(p1);
  }else{
    sumsig1 = dvec(1);
  }
  suma1 = dvec(p1);
  sumb1 = dvec(p1k);
  sumGamma1 = ivec(p1);
  sumAcceptGamma1 = ivec(p1);
  sumAcceptBeta1 = ivec(p1);

  if(ty2 == 1){
    sumsig2 = dvec(p2);
  }else{
    sumsig2 = dvec(1);
  }
  suma2 = dvec(p2);
  sumb2 = dvec(p2k);
  sumGamma2 = ivec(p2);
  sumAcceptGamma2 = ivec(p2);
  sumAcceptBeta2 = ivec(p2);

  if(ty3 == 1){
    sumsig3 = dvec(p3);
  }else{
    sumsig3 = dvec(1);
  }   
  suma3 = dvec(p3);
  sumb3 = dvec(p3k);
  sumGamma3 = ivec(p3);
  sumAcceptGamma3 = ivec(p3);
  sumAcceptBeta3 = ivec(p3);

  if(ty4 == 1){
    sumsig4 = dvec(p4);
  }else{
    sumsig4 = dvec(1);
  }   
  suma4 = dvec(p4);
  sumb4 = dvec(p4k);
  sumGamma4 = ivec(p4);
  sumAcceptGamma4 = ivec(p4);
  sumAcceptBeta4 = ivec(p4);

  if(ty5 == 1){
    sumsig5 = dvec(p5);
  }else{
    sumsig5 = dvec(1);
  }   
  suma5 = dvec(p5);
  sumb5 = dvec(p5k);
  sumGamma5 = ivec(p5);
  sumAcceptGamma5 = ivec(p5);
  sumAcceptBeta5 = ivec(p5);
  
  sumMeanZ = dvec(nk);
 /* Important to use GetRNGstate */
  GetRNGstate();
  /************* mcmc burnin ******************/
  for(i=0; i<(betaBurnin); i++){
     mcmcMix6d(meanZ,lastZ,nkZbd,sdev,
	       &ty0,&p0,a0,gb0,con0,cat0,sigma0,
	       &ty1,&p1,a1,gb1,con1,cat1,sigma1,
	       &ty2,&p2,a2,gb2,con2,cat2,sigma2,
	       &ty3,&p3,a3,gb3,con3,cat3,sigma3,
	       &ty4,&p4,a4,gb4,con4,cat4,sigma4,
	       &ty5,&p5,a5,gb5,con5,cat5,sigma5,acceptZ);

     standardize(meanZ,&n,&k);
     if(ndt > 0){
       switch(ty0){
       case 1 : bvsNormal(con0,meanZ,a0,b0,sigma0,gamma0,acceptGamma0,&pg[0],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p0,&k); break;
       case 2 : bvsBinom(cat0,meanZ,a0,b0,acceptBeta0,gamma0,acceptGamma0,&pg[0],invSigma0,beta0,&mhLambda,&n,&p0,&k); break;
       case 3 : bvsPoisson(cat0,meanZ,a0,b0,acceptBeta0,gamma0,acceptGamma0,&pg[0],invSigma0,beta0,&mhLambda,&n,&p0,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p0; j++){
	   gb0[ID] = b0[ID]*gamma0[j];
	   ID = ID + 1;
	 }
       } 
     }

     if(ndt > 1){
       switch(ty1){
       case 1 : bvsNormal(con1,meanZ,a1,b1,sigma1,gamma1,acceptGamma1,&pg[1],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p1,&k); break;
       case 2 : bvsBinom(cat1,meanZ,a1,b1,acceptBeta1,gamma1,acceptGamma1,&pg[1],invSigma0,beta0,&mhLambda,&n,&p1,&k); break;
       case 3 : bvsPoisson(cat1,meanZ,a1,b1,acceptBeta1,gamma1,acceptGamma1,&pg[1],invSigma0,beta0,&mhLambda,&n,&p1,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p1; j++){
	   gb1[ID] = b1[ID]*gamma1[j];
	   ID = ID + 1;
	 }
       }
     }

     if(ndt > 2){
       switch(ty2){
       case 1 : bvsNormal(con2,meanZ,a2,b2,sigma2,gamma2,acceptGamma2,&pg[2],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p2,&k); break;
       case 2 : bvsBinom(cat2,meanZ,a2,b2,acceptBeta2,gamma2,acceptGamma2,&pg[2],invSigma0,beta0,&mhLambda,&n,&p2,&k); break;
       case 3 : bvsPoisson(cat2,meanZ,a2,b2,acceptBeta2,gamma2,acceptGamma2,&pg[2],invSigma0,beta0,&mhLambda,&n,&p2,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p2; j++){
	   gb2[ID] = b2[ID]*gamma2[j];
	   ID = ID + 1;
	 }
       } 
     }

     if(ndt > 3){
       switch(ty3){
       case 1 : bvsNormal(con3,meanZ,a3,b3,sigma3,gamma3,acceptGamma3,&pg[3],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p3,&k); break;
       case 2 : bvsBinom(cat3,meanZ,a3,b3,acceptBeta3,gamma3,acceptGamma3,&pg[3],invSigma0,beta0,&mhLambda,&n,&p3,&k); break;
       case 3 : bvsPoisson(cat3,meanZ,a3,b3,acceptBeta3,gamma3,acceptGamma3,&pg[3],invSigma0,beta0,&mhLambda,&n,&p3,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p3; j++){
	   gb3[ID] = b3[ID]*gamma3[j];
	   ID = ID + 1;
	 }
       }        
     }
     
     if(ndt > 4){
       switch(ty4){
       case 1 : bvsNormal(con4,meanZ,a4,b4,sigma4,gamma4,acceptGamma4,&pg[4],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p4,&k); break;
       case 2 : bvsBinom(cat4,meanZ,a4,b4,acceptBeta4,gamma4,acceptGamma4,&pg[4],invSigma0,beta0,&mhLambda,&n,&p4,&k); break;
       case 3 : bvsPoisson(cat4,meanZ,a4,b4,acceptBeta4,gamma4,acceptGamma4,&pg[4],invSigma0,beta0,&mhLambda,&n,&p4,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p4; j++){
	   gb4[ID] = b4[ID]*gamma4[j];
	   ID = ID + 1;
	 }
       }        
     }

     if(ndt > 5){
       switch(ty5){
       case 1 : bvsNormal(con5,meanZ,a5,b5,sigma5,gamma5,acceptGamma5,&pg[5],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p5,&k); break;
       case 2 : bvsBinom(cat5,meanZ,a5,b5,acceptBeta5,gamma5,acceptGamma5,&pg[5],invSigma0,beta0,&mhLambda,&n,&p5,&k); break;
       case 3 : bvsPoisson(cat5,meanZ,a5,b5,acceptBeta5,gamma5,acceptGamma5,&pg[5],invSigma0,beta0,&mhLambda,&n,&p5,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p5; j++){
	   gb5[ID] = b5[ID]*gamma5[j];
	   ID = ID + 1;
	 }
       }        
     }     
  }

  /* acceptZ is added up during burnin; set acceptZ to zero */
  for(i=0; i<n; i++){
    acceptZ[i] = 0;
  }

  /* MCMC sampling */
  for(i=0; i<(betaDraw); i++){
    mcmcMix6d(meanZ,lastZ,nkZbd,sdev,
	       &ty0,&p0,a0,gb0,con0,cat0,sigma0,
	       &ty1,&p1,a1,gb1,con1,cat1,sigma1,
	       &ty2,&p2,a2,gb2,con2,cat2,sigma2,
	       &ty3,&p3,a3,gb3,con3,cat3,sigma3,
	       &ty4,&p4,a4,gb4,con4,cat4,sigma4,
	       &ty5,&p5,a5,gb5,con5,cat5,sigma5,acceptZ);    

     standardize(meanZ,&n,&k);
     remainder = (i % thin);
     if(remainder == 0){
       dvadd(sumMeanZ,meanZ,nk);
       nthin = nthin + 1;
     }

     if(ndt > 0){
       switch(ty0){
       case 1 : bvsNormal(con0,meanZ,a0,b0,sigma0,gamma0,acceptGamma0,&pg[0],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p0,&k); break;
       case 2 : bvsBinom(cat0,meanZ,a0,b0,acceptBeta0,gamma0,acceptGamma0,&pg[0],invSigma0,beta0,&mhLambda,&n,&p0,&k); break;
       case 3 : bvsPoisson(cat0,meanZ,a0,b0,acceptBeta0,gamma0,acceptGamma0,&pg[0],invSigma0,beta0,&mhLambda,&n,&p0,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p0; j++){
	   gb0[ID] = b0[ID]*gamma0[j];
	   ID = ID + 1;
	 }
       }     
       ivadd(sumGamma0,gamma0,p0);
       ivadd(sumAcceptGamma0,acceptGamma0,p0);
       ivadd(sumAcceptBeta0,acceptBeta0,p0);
       if(remainder == 0){
	 if(ty0 == 1){
	   dvadd(sumsig0,sigma0,p0); 
	 }
	 dvadd(suma0,a0,p0);
	 dvadd(sumb0,b0,p0k); 
       }
     }

     if(ndt > 1){
       switch(ty1){
       case 1 : bvsNormal(con1,meanZ,a1,b1,sigma1,gamma1,acceptGamma1,&pg[1],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p1,&k); break;
       case 2 : bvsBinom(cat1,meanZ,a1,b1,acceptBeta1,gamma1,acceptGamma1,&pg[1],invSigma0,beta0,&mhLambda,&n,&p1,&k); break;
       case 3 : bvsPoisson(cat1,meanZ,a1,b1,acceptBeta1,gamma1,acceptGamma1,&pg[1],invSigma0,beta0,&mhLambda,&n,&p1,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p1; j++){
	   gb1[ID] = b1[ID]*gamma1[j];
	   ID = ID + 1;
	 }
       }        
       ivadd(sumGamma1,gamma1,p1);
       ivadd(sumAcceptGamma1,acceptGamma1,p1);
       ivadd(sumAcceptBeta1,acceptBeta1,p1);
       if(remainder == 0){
	 if(ty1 == 1){
	   dvadd(sumsig1,sigma1,p1); 
	 }
	 dvadd(suma1,a1,p1);
	 dvadd(sumb1,b1,p1k); 
       }	
     }

     if(ndt > 2){
       switch(ty2){
       case 1 : bvsNormal(con2,meanZ,a2,b2,sigma2,gamma2,acceptGamma2,&pg[2],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p2,&k); break;
       case 2 : bvsBinom(cat2,meanZ,a2,b2,acceptBeta2,gamma2,acceptGamma2,&pg[2],invSigma0,beta0,&mhLambda,&n,&p2,&k); break;
       case 3 : bvsPoisson(cat2,meanZ,a2,b2,acceptBeta2,gamma2,acceptGamma2,&pg[2],invSigma0,beta0,&mhLambda,&n,&p2,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p2; j++){
	   gb2[ID] = b2[ID]*gamma2[j];
	   ID = ID + 1;
	 }
       }        
       ivadd(sumGamma2,gamma2,p2);
       ivadd(sumAcceptGamma2,acceptGamma2,p2);
       ivadd(sumAcceptBeta2,acceptBeta2,p2);
       if(remainder == 0){
	 if(ty2 == 1){
	   dvadd(sumsig2,sigma2,p2); 
	 }
	 dvadd(suma2,a2,p2);
	 dvadd(sumb2,b2,p2k); 
       }
     }

     if(ndt > 3){
       switch(ty3){
       case 1 : bvsNormal(con3,meanZ,a3,b3,sigma3,gamma3,acceptGamma3,&pg[3],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p3,&k); break;
       case 2 : bvsBinom(cat3,meanZ,a3,b3,acceptBeta3,gamma3,acceptGamma3,&pg[3],invSigma0,beta0,&mhLambda,&n,&p3,&k); break;
       case 3 : bvsPoisson(cat3,meanZ,a3,b3,acceptBeta3,gamma3,acceptGamma3,&pg[3],invSigma0,beta0,&mhLambda,&n,&p3,&k); break;
       case 4 : break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p3; j++){
	   gb3[ID] = b3[ID]*gamma3[j];
	   ID = ID + 1;
	 }
       } 
       ivadd(sumGamma3,gamma3,p3);
       ivadd(sumAcceptGamma3,acceptGamma3,p3);
       ivadd(sumAcceptBeta3,acceptBeta3,p3);
       if(remainder == 0){
	 if(ty3 == 1){
	   dvadd(sumsig3,sigma3,p3); 
	 }
	 dvadd(suma3,a3,p3);
	 dvadd(sumb3,b3,p3k); 
       }
     }

     if(ndt > 4){
       switch(ty4){
       case 1 : bvsNormal(con4,meanZ,a4,b4,sigma4,gamma4,acceptGamma4,&pg[4],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p4,&k); break;
       case 2 : bvsBinom(cat4,meanZ,a4,b4,acceptBeta4,gamma4,acceptGamma4,&pg[4],invSigma0,beta0,&mhLambda,&n,&p4,&k); break;
       case 3 : bvsPoisson(cat4,meanZ,a4,b4,acceptBeta4,gamma4,acceptGamma4,&pg[4],invSigma0,beta0,&mhLambda,&n,&p4,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p4; j++){
	   gb4[ID] = b4[ID]*gamma4[j];
	   ID = ID + 1;
	 }
       } 
       ivadd(sumGamma4,gamma4,p4);
       ivadd(sumAcceptGamma4,acceptGamma4,p4);
       ivadd(sumAcceptBeta4,acceptBeta4,p4);
       if(remainder == 0){
	 if(ty4 == 1){
	   dvadd(sumsig4,sigma4,p4); 
	 }
	 dvadd(suma4,a4,p4);
	 dvadd(sumb4,b4,p4k); 
       }
     }

     if(ndt > 5){
       switch(ty5){
       case 1 : bvsNormal(con5,meanZ,a5,b5,sigma5,gamma5,acceptGamma5,&pg[5],invSigma0,invSigmaBeta0,&invga0,&invgb0,&n,&p5,&k); break;
       case 2 : bvsBinom(cat5,meanZ,a5,b5,acceptBeta5,gamma5,acceptGamma5,&pg[5],invSigma0,beta0,&mhLambda,&n,&p5,&k); break;
       case 3 : bvsPoisson(cat5,meanZ,a5,b5,acceptBeta5,gamma5,acceptGamma5,&pg[5],invSigma0,beta0,&mhLambda,&n,&p5,&k); break;
       }
       ID = 0;
       for(h=0; h<k; h++){
	 for(j=0; j<p5; j++){
	   gb5[ID] = b5[ID]*gamma5[j];
	   ID = ID + 1;
	 }
       } 
       ivadd(sumGamma5,gamma5,p5);
       ivadd(sumAcceptGamma5,acceptGamma5,p5);
       ivadd(sumAcceptBeta5,acceptBeta5,p5);
       if(remainder == 0){
	 if(ty5 == 1){
	   dvadd(sumsig5,sigma5,p5); 
	 }
	 dvadd(suma5,a5,p5);
	 dvadd(sumb5,b5,p5k); 
       }
     }
  }
  PutRNGstate();
  
  scale = 1.0/nthin;
  dvscale(sumMeanZ,nk,scale);
  dvcopy(meanZ,sumMeanZ,nk);
  Free(sumMeanZ);

  if(ndt > 0){
    ivcopy(gamma0,sumGamma0,p0);
    ivcopy(acceptGamma0,sumAcceptGamma0,p0);
    ivcopy(acceptBeta0,sumAcceptBeta0,p0);
    dvscale(suma0,p0,scale);
    dvscale(sumb0,p0k,scale);
    dvcopy(a0,suma0,p0);
    dvcopy(b0,sumb0,p0k);
    if(ty0 == 1){
      dvscale(sumsig0,p0,scale);
      dvcopy(sigma0,sumsig0,p0); 
    }
    Free(gb0);
  }

  if(ndt > 1){
    ivcopy(gamma1,sumGamma1,p1);
    ivcopy(acceptGamma1,sumAcceptGamma1,p1);
    ivcopy(acceptBeta1,sumAcceptBeta1,p1);
    dvscale(suma1,p1,scale); 
    dvscale(sumb1,p1k,scale); 
    dvcopy(a1,suma1,p1);
    dvcopy(b1,sumb1,p1k);
    if(ty1 == 1){
      dvscale(sumsig1,p1,scale);
      dvcopy(sigma1,sumsig1,p1);
    }
    Free(gb1);
  }

  if(ndt > 2){
    ivcopy(gamma2,sumGamma2,p2);
    ivcopy(acceptGamma2,sumAcceptGamma2,p2);
    ivcopy(acceptBeta2,sumAcceptBeta2,p2);
    dvscale(suma2,p2,scale);
    dvscale(sumb2,p2k,scale);
    dvcopy(a2,suma2,p2);
    dvcopy(b2,sumb2,p2k);
    if(ty2 == 1){
      dvscale(sumsig2,p2,scale);
      dvcopy(sigma2,sumsig2,p2);
    }
    Free(gb2);
  }

  if(ndt > 3){
    ivcopy(gamma3,sumGamma3,p3);
    ivcopy(acceptGamma3,sumAcceptGamma3,p3);
    ivcopy(acceptBeta3,sumAcceptBeta3,p3);
    dvscale(suma3,p3,scale);
    dvscale(sumb3,p3k,scale);
    dvcopy(a3,suma3,p3);
    dvcopy(b3,sumb3,p3k);
    if(ty3 == 1){
      dvscale(sumsig3,p3,scale);
      dvcopy(sigma3,sumsig3,p3);
    }
    Free(gb3);
  }

  if(ndt > 4){
    ivcopy(gamma4,sumGamma4,p4);
    ivcopy(acceptGamma4,sumAcceptGamma4,p4);
    ivcopy(acceptBeta4,sumAcceptBeta4,p4);
    dvscale(suma4,p4,scale);
    dvscale(sumb4,p4k,scale);
    dvcopy(a4,suma4,p4);
    dvcopy(b4,sumb4,p4k);
    if(ty4 == 1){
      dvscale(sumsig4,p4,scale);
      dvcopy(sigma4,sumsig4,p4);
    }
    Free(gb4);
  }

  if(ndt > 5){
    ivcopy(gamma5,sumGamma5,p5);
    ivcopy(acceptGamma5,sumAcceptGamma5,p5);
    ivcopy(acceptBeta5,sumAcceptBeta5,p5);
    dvscale(suma5,p5,scale);
    dvscale(sumb5,p5k,scale);
    dvcopy(a5,suma5,p5);
    dvcopy(b5,sumb5,p5k);
    if(ty5 == 1){
      dvscale(sumsig5,p5,scale);
      dvcopy(sigma5,sumsig5,p5);
    }
    Free(gb5);
  }
  
  Free(sumGamma0);
  Free(sumAcceptGamma0);
  Free(sumAcceptBeta0);
  Free(sumGamma1);
  Free(sumAcceptGamma1);
  Free(sumAcceptBeta1);
  Free(sumGamma2);
  Free(sumAcceptGamma2);
  Free(sumAcceptBeta2); 
  Free(sumGamma3);
  Free(sumAcceptGamma3);
  Free(sumAcceptBeta3);
  Free(sumGamma4);
  Free(sumAcceptGamma4);
  Free(sumAcceptBeta4);
  Free(sumGamma5);
  Free(sumAcceptGamma5);
  Free(sumAcceptBeta5);   
  Free(suma0);
  Free(suma1);
  Free(suma2);
  Free(suma3);
  Free(suma4);
  Free(suma5);
  Free(sumb0);
  Free(sumb1);
  Free(sumb2);
  Free(sumb3);
  Free(sumb4);
  Free(sumb5);
  Free(sumsig0);
  Free(sumsig1);
  Free(sumsig2);
  Free(sumsig3);
  Free(sumsig4);
  Free(sumsig5);  
}
