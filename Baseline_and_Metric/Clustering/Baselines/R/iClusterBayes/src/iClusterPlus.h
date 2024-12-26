/* Qianxing Mo, Dan L. Duncan Cancer Center, Baylor College of Medicine */
/* Code programs for iCluster and giCluster
   2nd last updated 8/6/2011 
   last updated 12/7/2012, change logp to logpnull, logq to logqnull
*/
/*iCluster utility function */
/* utility.c derived from util.c; dataType and fillData() has been changed  */
/* 07/06/2011 */

typedef struct{
  int *type;
  int *p;  /* p variable */
  int *k;  /* k latent class of subject */
  int *c;  /* c class of categoriable variable (subject)  */
  double *alpha;
  double *beta;
  double *con; /* continous variable, should be a px1 vector */
  int *cat;   /* categorical or count variable, should be a px1 vector  */
  int *class;  /* the actual values of multinomial values, from 0 to (max-number-of-class - 1) */
  int *nclass; /* Number of class for each p variable, across samples */
  double *sigma2; /* for linear regression case */
}dataType;

int  *ivec(int len);
double *dvec(int len);

/*a <- b */
void dvcopy(double *a, double *b, int row);

/*a <- b */
void ivcopy(int *a, int *b, int row);

/* x = x - y*/
void dvsub(double * x, double *y, int size);

/* x = x + y */
void dvadd(double * x, double *y, int size);

/* x = x + y */
void ivadd(int * x, int *y, int size);

/* x = x*alpha */
void dvscale(double * x, int size, double alpha);

/*  inner product of x & y  */
double dvdot(double *x, double *y, int size);

/* x = alpha/x */
void dvinv(double * x, int size, double alpha);

/*A <- B[rowStart:rowEnd,colStart:colEnd] */
void dvsect(double *A, double *B, int Start, int End);

/* sum(vec^2) */
double vecsum2(double *rvec, int n);

/* B = inverse(A) */
/* Note A is also change (see dgesv manual) */
/* should be very careful using this function */
void invsqm(double *B, double *A, int *n);

/* B = inverse(A) */
/* Note A is NOT change */
/* invsqm2 is theoretically slowed than invsqm */
void invsqm2(double *B, double *A, int *n);

/* x = x + y */
void dmadd(double ** x, double **y, int row, int col);

/* x = x - y */
void dmsub(double ** x, double **y, int row, int col);

/* x = scale*x */
void dmscale(double ** x, int row, int col, double scale);

/* x = apply(y,1,sum) */
void dmrowsum(double * x, double **y, int row, int col);

/*A <- B[rowStart:rowEnd,colStart:colEnd] */
void dmsect(double *A, double *B, int rowB,int rowStart, int rowEnd, int colStart, int colEnd);

/*B[rowStart:rowEnd,colStart:colEnd] = A */
void dmreplace(double *B, double *A, int rowB,int rowStart, int rowEnd, int colStart, int colEnd);

/*A <- B */
void dmcopy(double **A, double **B, int row, int col);

/* A = b, A is a row matrix */
void dvtom(double **A, double *b, int brow, int bcol);

/* A = b, A is a row matrix, b is a vector */
void ivtom(int **A, int *b, int brow, int bcol);

/*b = A, b is column vector */
void dmtov( double *b,double **A, int arow, int acol);

/*create a diag matrix with a single value, diag(m)=val*/
void diagm(double *m, int row,int val);

/*create a diag matrix using a vector, m is row by row  */
/* m = diag(val) */
void diagmv(double *m, int row,double *val);

/*get the diag elements of a square matrix */
void diagv(double *v, double *m, int row);

/*diag elements of m + scale */
void diagplus(double *m, int nrow, double scale);

/*diag elements of m + scale vector */
void diagplusv(double *m, int nrow, double *scale);

/*row matrix */
double **drowm(int row, int col);

/*integer row matrix */
int **irowm(int row, int col);

/*column matrix */
double **dcolm(int row, int col);

void dmfree(double **v,int row);

void difree(int **v,int row);

void imfree(int **v,int row);

/*if abs(m[i][j]) < val, m[i][k] = val  */
void editm(double **m, int row, int col, double val);

void printmatrix(double **m, int row,int col);

void printvec(double *v, int len);

/* v1= 1/abs(v2) */
void  fabsinv(double *v1,double *v2, int len);

/* v = t(m)*, v is column matrix */
void dmtranv(double *v, double **m, int row, int col);

/* m1 = t(m2), m1 is m2col by m2row; m2 is m2row by m2col */
void dmtranm(double **m1, double **m2, int m2row, int m2col);

void fillData(dataType *dt,double *alpha,double *beta,double *con,int *cat,int *class,int *nclass,
	      double *sigma2,int *p,int *k,int *c,int *type);

/* used for iClusterBayes, which only modle normal, poisson and binomial data */
void fillData3(dataType *dt,double *alpha,double *beta,double *con,int *cat,
	       double *sigma2,int *p,int *k,int *type);

/*iCluster C programs -- implement Lasso, Enet, Group Lasso, Fussed Lasso, Group Fussed Lasso. 
  Last updated: 10/11/2010
  Qianxing Mo (qmo@bcm.edu), Dan L. Duncan Cancer Center, Baylor College of Medicine
*/

/*X is the t(X) of the original R code */
void  Mstep_glasso(int *p, int *k, int *n,double *B, double *X, double *Phivec, 
		   double *EZZt, double *EZ, double *lam, double *eps2);

/*X is the t(X) of the original R code */
void  Mstep_lasso(int *p, int *k, int *n,double *B, double *X, double *Phivec, 
		  double *EZZt, double *EZ, double *lam, double *eps2);

/*X is the t(X) of the original R code */
void  Mstep_enet(int *p, int *k, int *n,double *B, double *X, double *Phivec, 
		 double *EZZt, double *EZ, double *lam, double *eps2);

void  Mstep_flasso(int *p, int *k,double *B, double *Phivec, double *EXZt,
		   double *EZZt, double *lam, double *eps2, int *id,int *idlen);

/* w = dvec(n);  eigen values
   z= dvec(n*n); eigen vectors
*/
void eigen(double *A, int *row, double *w, double *z);

void lyap(double *B, double *P, double *Q, double *C, int *m, int *n);

void  Mstep_gflasso(int *p,int *k,double *B,double *Phivec,double *EXZt,
		    double *EZZt,double *lam,double *eps2,int *id,int *idlen);

/* dim(B) = p,k; dim(X) = n,p */
void iClusterCore(int *p, int *k, int *n, double *xtxdiag, double *X,double *B,double *EZ,double *EZZt,
		  double *Phivec, double *dif, int *iter, int *pvec,
		  double *lambda,double *eps,double *eps2, int *maxiter, int *lenT,int *method, int *ID, 
		  int *lenID);

/*iClusterPlus.c is the merged file of glmnetc.c and mixmcmc.c */
/* glmnet C wrapper functions; continous work on elnet.c; handle mixed type of categorical data by modifying lognetC and lognetBath */
/* Qianxing Mo, Dan L. Duncan Cancer Center, Baylor College of Medicine */
/* Last update: 7/13/2011  */

/* Fortran function of glmnet */
void elnet_(int *ka,double *alpha,int *nobs,int *nvars,double *x,double *y,double *weights,
	    int *jd,double *vp,int *ne,int *nx,int *nlam,double *flmin,double *ulam,double *thresh,
	 int *isd,int *maxit,int *lmu,double *a0,double *ca,int *ia,int *nin,double *rsq,double *alm,
	    int *nlp,int *jerr);

void fishnet_(double *parm,int *no,int *ni,double *x,double *y,double *g,double *w,int *jd,double *vp,
	      int *ne,int *nx,int *nlam,double *flmin,double *ulam,double *thr,int *isd,
	      int *maxit,int *lmu,double *a0,double *ca,int *ia,int *nin,double *dev0,double *dev,
	      double *alm,int *nlp,int *jerr);

void lognet_(double *parm,int *no,int *ni,int *nc,double *x,double *y,double *g,int *jd,
	     double *vp,int *ne,int *nx,int *nlam,double *flmin,double *ulam,double *thr,int *isd,
	     int *maxit,int *kopt,int *lmu,double *a0,double *ca,int *ia,int *nin,double *dev0,double *dev,
	      double *alm,int *nlp,int *jerr);

void getbeta(double *beta,int *df,int *nin,int *nvars,int *ia,double *ca);
void getbetaMult(double *beta,int *df,int *nin,int *nvars,int *nc,int *ia,double *ca);

/* elnetC call glmnet Fortran function elnet_; it work only when lambda is a scalar (nlam==1)  */
/* nobs = dim(x)[1], nvar=dim(x)[2], ulam=lambda */
/* vp = penalty.factor = rep(1,nvars),weights=rep(1,nobs) */
/* a0, beta[nvars],df,rsq are output; rsq is the dev.ratio = 1 - dev/nulldev */
/* rsq = dev.ratio = (1-dev/nulldev) */
void elnetC(double *a0,double *beta,int *df,double *x,double *y,int *nobs,int *nvars,
	    double *alpha,double *lambda,double *rsq);

/* a0[p],beta[p][k],df[p],x[nobs][k],y[nobs][p], y = a0 + beta'x */
void elnetBatch(double *a0,double *beta,double *sigma2,double *x,double *y,int *nobs,int *k,int *p,
		double *alpha,double *lambda);

/* used to calculate sum of deviance and null deviance */
void elnetBatchDev(double *a0,double *beta,double *sigma2,double *x,double *y,int *nobs,int *k,int *p,
		double *alpha,double *lambda,double *sumdev0,double *sumdev);


/* fishnetC call glmnet Fortran function elnet_; it work only when lambda is a scalar (nlam==1)  */
/* nobs = dim(x)[1], nvar=dim(x)[2], ulam=lambda */
/* vp = penalty.factor = rep(1,nvars),weights=rep(1,nobs) */
/* a0, beta[nvars],df,dev0 (null deviance) and dev (dev.ratio) are output;  */
/* dev0 is null deviance; dev = dev.ratio = (1-dev/nulldev) */
void fishnetC(double *a0,double *beta,int *df,double *x,double *y,int *nobs,int *nvars,
	      double *alpha,double *lambda,double *dev0,double *dev);

void fishnetBatch(double *a0,double *beta,double *x,double *y,int *nobs,int *k,int *p,
		  double *alpha,double *lambda);

/* used to calculate sum of deviance and null deviance */
void fishnetBatchDev(double *a0,double *beta,double *x,double *y,int *nobs,int *k,int *p,
		  double *alpha,double *lambda,double *sumdev,double *sumnulldev);

/* lognetC call glmnet Fortran function lognet_; it work only when lambda is a scalar (nlam==1)  */
/* y should be factor starting from 0; perform as.numeric(as.factor(y)) - 1 in R*/
/* nobs = dim(x)[1], nvar=dim(x)[2], ulam=lambda */
/* vp = penalty.factor = rep(1,nvars),weights=rep(1,nobs) */
/* a0, beta[nvars] and df are output;  */
/* if family = binomial; a0[1], beta[nvars]; if family = multinomial, a0[nclass],beta[nclass*nvars] */
/* nclass is a scalar */
/* dev0 is null deviance; dev = dev.ratio = (1-dev/nulldev) */
void lognetC(double *a0,double *beta,int *df,double *x,int *y,int *nobs,int *nvars,
	     double *alpha,double *lambda,int *nclass,int *family,double *dev0,double *dev);

/* a0[p],beta[p][k],df[p],x[nobs][k],y[nobs][p], y = a0 + beta'x */
/* nclass is vector, indicating the number of class for the columns of y */
/* maxclass is the maximum of nclass */
void lognetBatch(double *a0,double *beta,double *x,int *y,int *nobs,int *k,int *p,
		 double *alpha,double *lambda,int *nclass,int *maxclass,int *family);

/* used to calculate sum of deviance and null deviance */
/* a0[p],beta[p][k],df[p],x[nobs][k],y[nobs][p], y = a0 + beta'x */
/* nclass is vector, indicating the number of class for the columns of y */
/* maxclass is the maximum of nclass */
void lognetBatchDev(double *a0,double *beta,double *x,int *y,int *nobs,int *k,int *p,
		 double *alpha,double *lambda,int *nclass,int *maxclass,int *family,
		 double *dev0,double *sumdev);

/* this function only work when lmu == 1 (nlam==1) */
/* get the coefficients from the elnet fortran program */
/* nin == maxnin when lambda is a scalar; number of compressed coefs for each solution */
/* beta is a vector with nvars elements; nx is not used here  */
void getbeta(double *beta,int *df,int *nin,int *nvars,int *ia,double *ca);

/* this function only work when lmu == 1 (nlam==1) */
/* get the coefficients from the lognet fortran program */
/* nin == maxnin when lambda is a scalar; number of compressed coefs for each solution */
/* beta is a matrix(nvars x nc); nx is not used here  */
void getbetaMult(double *beta,int *df,int *nin,int *nvars,int *nc,int *ia,double *ca);

double nulldev(double *y,int n);

/* the following code is from mixmcmc.c */
/* continous work on mixpca.c; handle mixed type of categorical data by modifying logMult,metroMix and mcmcXmix */ 
/* Qianxing Mo, Dan L. Duncan Cancer Center */
/*Last update: 7/5/2011  */

/*  logpnull <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
    logqnull <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta))) */
double logpnull(double etai);

double logqnull(double etai);

/* log Binomial likelihood for subject i */
/* alpha is p */
/* beta is p by k */
/* x is 1 by p  */
/* zi is 1 by k  */
double logBinom(double *zi,double *alpha,double *beta,int *xi,int *p,int *k);

/* log Binomial likelihood for subject i to n */
/* alpha is p */
/* beta is p by k */
/* X is n by p  */
/* Z is n by k  */
void logBinomAll(double *LogLike,double *Z,double *alpha,double *beta,int *X,int *n, int *p,int *k);

/* log Poisson likelihood for subject i */
/* alpha is p */
/* beta is p by k */
/* x is 1 by p  */
/* zi is 1 by k  */
double logPoisson(double *zi,double *alpha,double *beta,int *xi,int *p,int *k);

/* log Poisson likelihood for subject i to n */
/* alpha is p */
/* beta is p by k */
/* X is n by p  */
/* Z is n by k  */
void logPoissonAll(double *LogLike,double *Z,double *alpha,double *beta,int *X,int *n,int *p,int *k);

/* log normal likelihood for subject i */
/* alpha is p */
/* beta is p by k */
/* yi is 1 by p  */
/* zi is 1 by k  */
/* sigma2 is p */
double logNorm(double *zi,double *alpha,double *beta,double *sigma2,double *yi,int *p,int *k);

/* log normal likelihood for subject i to n */
/* alpha is p */
/* beta is p by k */
/* Y is n by p  */
/* Z is n by k  */
/* sigma2 is p */

void logNormAll(double *LogLike,double *Z,double *alpha,double *beta,double *sigma2,double *Y,int *n,int *p,int *k);

int indicator(int x, int y);

/* log multinomial likelihood for subject i */
/*
xi = 1xp; for the data matrix X[n,p], each column of X[,j] must be transformed by 
doing as.numeric(as.factor(X[,j])) in R eviroment 
beta = p  x (k x C) is the format in R; beta must be transposed using t(beta) before passing to this function 
alpha = p x C; alpha doesn't need to be transposed. 
zi = 1 x k
class is the class ID, e.g., 0,1,2
nclass = 1 x p, is the number of class for each p[i], i = 0, ..., p-1
C is the maximum number of class in all p markers 
*/

double logMult(double *zi,double *alpha,double *beta,int *xi,int *class,int *nclass,int *p,int *C,int *K);

/* log multinomial likelihood for subject i to n*/
/*
X = n x p; for the data matrix X[n,p], each column of X[,j] must be transformed by 
doing as.numeric(as.factor(X[,j])) in R eviroment 
beta = p  x (k x C) is the format in R; beta must be transposed using t(beta) before passing to this function 
alpha = p x C; alpha doesn't need to be transposed. 
Z = n x k
class = 1 x C; is the class ID, e.g., 0,1,2
nclass = 1 x p; is the number of class for each p[i], i = 0, ..., p-1
C is the maximum number of class in all p markers 
*/
void logMultAll(double *LogLike,double *Z,double *alpha,double *beta,int *X,int *class,int *nclass,int *n,int *p,int *C,int *K);

void metroBinom(double *zi,double *alpha,double *beta,int *xi,int *p,int *k,double *sigma,double *newz);

void metroPoisson(double *zi,double *alpha,double *beta,int *xi,int *p,int *k,double *sigma,double *newz);

/* beta = p  x (k x C) is the format in R; beta must be transposed using t(beta) before 
   passing to this function */
void metroMult(double *zi,double *alpha,double *beta,int *xi,int *class,int *nclass,int *p,int *c,int *k,double *sigma,double *newz);

void metroNormal(double *zi,double *mu,double *beta, double *sigma2,double *yi,int *p,int *k,double *sigma,double *newz);

/*if tryz is accepted, accept = accept + 1 */
void metroMix(double *zi,dataType *dt1,dataType *dt2,dataType *dt3,dataType *dt4, double *sigma,
	      int *ndt,int *accept);

/* 11/14/2013 Q Mo rewrite matrix dynamic allocation and free code to avoid compiling warning on windows */
/* 2/28/2013, Q Mo fixed the potential memory releasing problem when data <= 2 */
/* if the data type is multinomial, beta = p  x (k x C) is the format in R; 
   beta must be transposed using t(beta) before passing to this function */
/* lastZ is changed after running the function */
/* meanZ is the mean of latent variable Z, lastZ is the last draw of Z;
   n: number of subjects; k: number of clusters; sdev - used to control random walk variation; 
   ndt: number of data types; typ; data type, 1 is normal, 2 is binary, 3 poisson, 4 multinomial; 
   c0-3: the max number of classes of the response if it is categorical;
   con0-3: continuous response; cat0-3: categorical response; 
   class0-3: the actual class of the categorical response; 
   nclass0-3: the number of class for categorical response
   accept: n x 1 vector (should be initialized as 0) to track if the proposal latent variable is accepted 
*/
void mcmcMix(double *meanZ,double *lastZ,int *n,int *k,int *burnin,int *draw,double *sdev,int *ndt,
	     int *ty0,int *p0,int *c0,double *a0,double *b0,double *con0,int *cat0,int *class0,int *nclass0,double *sigma0,
	     int *ty1,int *p1,int *c1,double *a1,double *b1,double *con1,int *cat1,int *class1,int *nclass1,double *sigma1,
	     int *ty2,int *p2,int *c2,double *a2,double *b2,double *con2,int *cat2,int *class2,int *nclass2,double *sigma2,
	     int *ty3,int *p3,int *c3,double *a3,double *b3,double *con3,int *cat3,int *class3,int *nclass3,double *sigma3,
	     int *accept);

/*
meanZ, lastZ are the output
initZ is the input;
*/
void mcmcBinom(double *meanZ,double *lastZ,double *alpha,double *beta,int *X, double *sigma, double *initz,
	       int *n,int *p,int *k,int *burnin,int *draw);
/*
meanZ, lastZ are the output
initZ is the input;
*/
void mcmcPoisson(double *meanZ,double *lastZ,double *alpha,double *beta,int *X, double *sigma, double *initz,
		 int *n,int *p,int *k,int *burnin,int *draw);

/*sigma2 is the estimated variances */
void mcmcNormal(double *meanZ,double *lastZ,double *mu,double *beta,double *sigma2,double *Y, double *sigma, 
		double *initz, int *n,int *p,int *k,int *burnin,int *draw);

/* if the data type is multinomial, beta = p  x (k x C) is the format in R; 
   beta must be transposed using t(beta) before passing to this function */
void mcmcMult(double *meanZ,double *lastZ,double *alpha,double *beta,int *X, int *class,int *nclass,
	      double *sigma,double *initz,int *n,int *p, int *c, int *k,int *burnin,int *draw);

/*scale(Z,center=TRUE,scale=TRUE) */
/* output Z will be standardized */
void standardize(double *Z,int *n,int *k);

/*
algorithm: 
covariance = inverse(precision)
cholesky decomposition of covariance: covariance = LL'
x = mean + Lz
elements in z are iid N(0,1)
*/

/* multivariate normal distribution with mean mu and covariance cov; vec is result. */
void rmvnormal(double *vec, double *mu, double *cov, int *n);

/* X is n x p  */
/* Z is n x k  */
/* alpha is p */
/* beta is p x k */
/* sigma2, gamma, accept_gamma are p x 1 */
/* invSigma0 is k+1 x k+1 */
/* invSigmaBeta0 is k+1 */
/* invga0 and invgb0 are scalar */
void bvsNormal(double *X,double *Z,double *alpha,double *beta,double *sigma2,int *gamma,int *accept_gamma,double *prior_ga1,double *invSigma0,double *invSigmaBeta0,double *invga0,double *invgb0,int *n,int *p,int *k);

/* joint sampling gamma and beta */
/* X is n x p  */
/* Z is n x k  */
/* alpha is p */
/* beta is p x k */
/* gamma, accept_gamma are p */
/* invSigma0 is k+1 x k+1 */
/* beta0 is k+1 */
/* lambda_mh is scalar */

void bvsPoisson(int *X,double *Z,double *alpha,double *beta,int *accept_beta,int *gamma,int *accept_gamma,double *prior_ga1,double *invSigma0,double *beta0,double *lambda_mh, int *n,int *p,int *k);

/* jointly sample alphaBeta and gamma */
/* X is n x p  */
/* Z is n x k  */
/* alpha is p */
/* beta is p x k */
/* gamma, accept_gamma are p */
/* invSigma0 is k+1 x k+1 */
/* beta0 is k+1 */
/* lambda_mh is scalar */
void bvsBinom(int *X,double *Z,double *alpha,double *beta,int *accept_beta,int *gamma,int *accept_gamma,double *prior_ga1,double *invSigma0,double *beta0,double *lambda_mh, int *n,int *p,int *k);

/*if tryz is accepted, accept = accept + 1; used for iClusterBayes */
void metroMix6d(double *zi,dataType *dt1,dataType *dt2,dataType *dt3,dataType *dt4,dataType *dt5,dataType *dt6,double *sigma,
		int *ndt,int *accept);

/* This function is modified from mcmcMix in the iClusterPlus.c; this fun can handle 6 data sets */
/* lastZ is changed after running the function */
/* meanZ is the mean of latent variable Z, lastZ is the last draw of Z; when ndraw=1,meanZ==lastZ,
   nkbd[0] is n; nkbd[1] is k; nkbd[2] is burnin; nkbd[3] is draw; nkbd[4] is ndt;
   n: number of subjects; k: number of clusters; sdev - used to control random walk variation; 
   ndt: number of data types; typ; data type, 1 is normal, 2 is binary, 3 poisson,(4 multinomial--not modeled in iClusterBayes); 
   a0-5 is p dimensional,b0-5[p][k],con0-5[n][p] or cat0-5[n][p], sigma0-5[p] 
   c0-3: the max number of classes of the response if it is categorical;
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
	       int *accept);

/* 11/14/2013 Q Mo rewrite matrix dynamic allocation and free code to avoid compiling warning on windows */
/* 2/28/2013, Q Mo fixed the potential memory releasing problem when data <= 2 */
/* if the data type is multinomial, beta = p  x (k x C) is the format in R; 
   beta must be transposed using t(beta) before passing to this function */
/* lastZ is changed after running the function */
/* meanZ is the mean of latent variable Z, lastZ is the last draw of Z;
   n: number of subjects; k: number of clusters; sdev - used to control random walk variation; 
   ndt: number of data types; typ; data type, 1 is normal, 2 is binary, 3 poisson, 4 multinomial; 
   c0-3: the max number of classes of the response if it is categorical;
   con0-3: continuous response; cat0-3: categorical response; 
   class0-3: the actual class of the categorical response; 
   nclass0-3: the number of class for categorical response
   accept: n x 1 vector (should be initialized as 0) to track if the proposal latent variable is accepted 
   nkZbd =c(n,k,Zburnin,Zdraw,ndt,thin)
   thin is used to control the dense of sumMeanZ and suma0-3,sumb0-3; if betaDraw % thin == 0, add meanZ,a0-3,b0-3,respectively.
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
	       int *gamma4,int *acceptBeta4,int *acceptGamma4,int *gamma5,int *acceptBeta5,int *acceptGamma5);
