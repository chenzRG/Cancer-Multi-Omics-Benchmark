/* Qianxing Mo,Department of Biostatistics & Bioinformatics, Moffitt Cancer Center */
/* Code programs for iCluster and giCluster
   2nd last updated 8/6/2011 
   last updated 12/7/2012, change logp to logpnull, logq to logqnull
*/
/*iCluster utility function */
/* utility.c derived from util.c; dataType and fillData() has been changed  */
/* 07/06/2011 */

#define USE_FC_LEN_T
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h> 
#include <R_ext/Lapack.h>
#include "iClusterPlus.h"

#ifndef FCONE
#define FCONE
#endif

/* #include <Rinternals.h>
   #include <R_ext/Rdynload.h>
   #include <R_ext/Utils.h>
*/


#define dim(A) INTEGER(coerceVector(getAttrib(A,R_DimSymbol),INTSXP))

int  *ivec(int len){
  int *v;
  v = (int *)Calloc(len, int);
  if(v == NULL){
    error("Error: fail to allocate memory space.\n");
  }
  return(v);
}

double *dvec(int len){
  double *v;
  v = (double *)Calloc(len, double);
  if(v == NULL){
    error("Error: fail to allocate memory space.\n");
  }
  return(v);
}

/*a <- b */
void dvcopy(double *a, double *b, int row){
  int i;
  for(i=0; i<row; i++){
    a[i] = b[i];
  }
}


/*a <- b */
void ivcopy(int *a, int *b, int row){
  int i;
  for(i=0; i<row; i++){
    a[i] = b[i];
  }
}

/* x = x - y*/
void dvsub(double * x, double *y, int size) {
  int i;
  for(i = 0; i < size; i++) {
    x[i] = x[i] - y[i];
  }
}
/* x = x + y */
void dvadd(double * x, double *y, int size) {
  int i;
  for(i = 0; i < size; i++) {
    x[i] = x[i] + y[i];
  }
}
/* x = x*alpha */
void dvscale(double * x, int size, double alpha) {
  int i;
  for(i = 0; i < size; i++) {
    x[i] = alpha*x[i];
  }
}

/*  inner product of x & y  */
double dvdot(double *x, double *y, int size) {
  int i;
  double sum;
  sum = 0;
  for(i = 0; i < size; i++) {
    sum = sum + x[i] * y[i];
  }
  return(sum);
}

/* x = alpha/x */
void dvinv(double * x, int size, double alpha) {
  int i;
  for(i = 0; i < size; i++) {
    x[i] = alpha/x[i];
  }
}

/*A <- B[rowStart:rowEnd,colStart:colEnd] */
void dvsect(double *A, double *B, int Start, int End){
  int i,id;
  id = 0;
  for(i = Start; i <= End; i++){
    A[id] = B[i];
    id++;
  }
}

/* sum(vec^2) */
double vecsum2(double *rvec, int n){
  int i;
  double value=0;
  for(i=0; i< n; i++){
    value += rvec[i]*rvec[i];
  }
  return value;
}

/* B = inverse(A) */
/* Note A is also change (see dgesv manual) */
/* should be very careful using this function */
void invsqm(double *B, double *A, int *n){
  int *IPIV;
  int i, j,id,INFO;
  IPIV = (int *)Calloc(*n, int);
  for(i=0; i<(*n); i++){
    for(j=0; j<(*n); j++){
      id = i*(*n)+j;
      if(i != j){
	B[id] = 0;
      }else{
	B[id] = 1;
      }
    }
  }
  F77_CALL(dgesv)(n,n,A,n,IPIV,B,n,&INFO);
  Free(IPIV);
}

/* B = inverse(A) */
/* Note A is NOT change */
/* invsqm2 is theoretically slowed than invsqm */
void invsqm2(double *B, double *A, int *n){
  int *IPIV;
  double *tempA;
  int i, j,id,nn,INFO;
  nn = (*n)*(*n);
  tempA = dvec(nn);
  dvcopy(tempA,A,nn);
  IPIV = (int *)Calloc(*n, int);
  for(i=0; i<(*n); i++){
    for(j=0; j<(*n); j++){
      id = i*(*n)+j;
      if(i != j){
	B[id] = 0;
      }else{
	B[id] = 1;
      }
    }
  }
  F77_CALL(dgesv)(n,n,tempA,n,IPIV,B,n,&INFO);
  Free(IPIV);
  Free(tempA);
}

/* x = x + y */
void dmadd(double ** x, double **y, int row, int col) {
  int i, j;
  for(i = 0; i < row; i++) {
    for(j=0; j < col; j++){
      x[i][j] = x[i][j] + y[i][j];
    }
  }
}

/* x = x - y */
void dmsub(double ** x, double **y, int row, int col) {
  int i, j;
  for(i = 0; i < row; i++) {
    for(j=0; j < col; j++){
      x[i][j] = x[i][j] - y[i][j];
    }
  }
}

/* x = scale*x */
void dmscale(double ** x, int row, int col, double scale) {
  int i, j;
  for(i = 0; i < row; i++) {
    for(j=0; j < col; j++){
      x[i][j] = scale * x[i][j];
    }
  }
}

/* x = apply(y,1,sum) */
void dmrowsum(double * x, double **y, int row, int col){
  int i, j;
  for(i = 0; i < row; i++) {
    x[i] = 0;
    for(j=0; j < col; j++){
      x[i] = x[i] + y[i][j];
    }
  }
}


/*A <- B[rowStart:rowEnd,colStart:colEnd] */
void dmsect(double *A, double *B, int rowB,int rowStart, int rowEnd, int colStart, int colEnd){
  int i, j,id;
  id = 0;
  for(i = (colStart); i <=(colEnd); i++){
    for(j = (rowStart); j <= (rowEnd); j++){
      A[id] = B[i*(rowB)+j];
      id++;
    }
  }
}


/*B[rowStart:rowEnd,colStart:colEnd] = A */
void dmreplace(double *B, double *A, int rowB,int rowStart, int rowEnd, int colStart, int colEnd){
  int i, j,id;
  id = 0;
  for(i = (colStart); i <= (colEnd); i++){
    for(j = (rowStart); j <= (rowEnd); j++){
      B[i*(rowB)+j] = A[id];
      id++;
    }
  }
}


/*A <- B */
void dmcopy(double **A, double **B, int row, int col){
  int i, j;
  for(i=0; i<row; i++){
    for(j=0; j<col; j++){
      A[i][j] = B[i][j];
    }
  }
}

/* A = b, A is a row matrix */
void dvtom(double **A, double *b, int brow, int bcol){
  int i, j;
  for(i=0; i<bcol; i++){
    for(j=0; j<brow;j++){
      A[j][i] = b[i*brow+j];
    }
  }
}

/* A = b, A is a row matrix, b is a vector */
void ivtom(int **A, int *b, int brow, int bcol){
  int i, j;
  for(i=0; i<bcol; i++){
    for(j=0; j<brow;j++){
      A[j][i] = b[i*brow+j];
    }
  }
}

/*b = A, b is column vector */
void dmtov( double *b,double **A, int arow, int acol){
  int i, j;
  for(i=0; i<acol; i++){
    for(j=0; j<arow;j++){
      b[i*arow+j] = A[j][i];
    }
  }
}

/*create a diag matrix with a single value, diag(m)=val*/
void diagm(double *m, int row,int val){
  int i, j,id;
  for(i=0; i<row; i++){
    for(j=0; j<row; j++){
      id = i*row+j;
      if(i != j){
	m[id] = 0;
      }else{
	m[id] = val;
      }
    }
  }
}

/*create a diag matrix using a vector, m is row by row  */
/* m = diag(val) */
void diagmv(double *m, int row,double *val){
  int i, j,id;
  for(i=0; i<row; i++){
    for(j=0; j<row; j++){
      id = i*row+j;
      if(i != j){
	m[id] = 0;
      }else{
	m[id] = val[i];
      }
    }
  }
}


/*get the diag elements of a square matrix */
void diagv(double *v, double *m, int row){
  int i, j,id;
  for(i=0; i<row; i++){
    for(j=0; j<row; j++){
      id = i*row+j;
      if(i == j){
	v[i] = m[id];
      }
    }
  }
}


/*diag elements of m + scale */
void diagplus(double *m, int nrow, double scale){
  int i;
  for(i=0; i<nrow;i++){
    m[i+nrow*i] += scale;
  }
}

/*diag elements of m + scale vector */
void diagplusv(double *m, int nrow, double *scale){
  int i;
  for(i=0; i<nrow;i++){
    m[i+nrow*i] += scale[i];
  }
}

/*row matrix */
double **drowm(int row, int col){
  int i;
  double **m;
  m = (double **)Calloc(row, double *);
  if(m == NULL){
    error("Error: fail to allocate memory space.\n");
  }
  for(i=0; i<row; i++) {
    m[i] = (double *)Calloc(col, double);
    if(m[i] == NULL){
      error("Error: fail to allocate memory space.\n");
    }
  }
  return(m);
}

/*integer row matrix */
int **irowm(int row, int col){
  int i;
  int **m;
  m = (int **)Calloc(row, int *);
  if(m == NULL){
    error("Error: fail to allocate memory space.\n");
  }
  for(i=0; i<row; i++) {
    m[i] = (int *)Calloc(col, int);
    if(m[i] == NULL){
      error("Error: fail to allocate memory space.\n");
    }
  }
  return(m);
}

/*column matrix */
double **dcolm(int row, int col){
  int i;
  double **m;
  m = (double **)Calloc(col, double *);
  if(m == NULL){
    error("Error: fail to allocate memory space.\n");
  }
  for(i=0; i<col; i++) {
    m[i] = (double *)Calloc(row, double);
    if(m[i] == NULL){
      error("Error: fail to allocate memory space.\n");
    }
  }
  return(m);
}

void dmfree(double **v,int row){
  int i;
  for(i=0; i<row; i++){
    Free(v[i]);
  }
  Free(v);  
}

void difree(int **v,int row){
  int i;
  for(i=0; i<row; i++){
    Free(v[i]);
  }
  Free(v);  
}

void imfree(int **v,int row){
  int i;
  for(i=0; i<row; i++){
    Free(v[i]);
  }
  Free(v);
}

/*if abs(m[i][j]) < val, m[i][k] = val  */
void editm(double **m, int row, int col, double val){
  int i,j;
  for(i=0; i<row; i++){
    for(j=0; j<col; j++){
      if(fabs(m[i][j]) < val){
	m[i][j] = val;
      }
    }
  }  
}

void printmatrix(double **m, int row,int col){
  int i,j;
  for(i=0; i<row; i++){
    for(j=0; j<col; j++){
      Rprintf("%f ",m[i][j]);
    }
    Rprintf("\n");
  }
}

void printvec(double *v, int len){
  int i;
  for(i=0; i<len; i++){
    Rprintf("%f ",v[i]);
  }
  Rprintf("\n\n");
}

/* v1= 1/abs(v2) */
void  fabsinv(double *v1,double *v2, int len){
  int i;
  for(i=0; i<len; i++){
    v1[i] = 1.0/fabs(v2[i]);
  }
}

/* v = t(m)*, v is column matrix */
void dmtranv(double *v, double **m, int row, int col){
  int i, j,id=0;
  for(i=0; i<row; i++){
    for(j=0; j<col; j++){
      v[id] = m[i][j];
      id = id + 1;
    }
  }
}

/* m1 = t(m2), m1 is m2col by m2row; m2 is m2row by m2col */
void dmtranm(double **m1, double **m2, int m2row, int m2col){
  int i, j;
  for(i=0; i<m2row; i++){
    for(j=0; j<m2col; j++){
      m1[j][i] = m2[i][j];
    }
  }
}

void fillData(dataType *dt,double *alpha,double *beta,double *con,int *cat,int *class,int *nclass,
	      double *sigma2,int *p,int *k,int *c,int *type){
  dt->type = type;
  dt->p = p;
  dt->k = k; 
  dt->c = c;
  dt->alpha = alpha;
  dt->beta = beta;
  dt->con = con;
  dt->cat = cat;
  dt->class = class;
  dt->nclass = nclass;
  dt->sigma2 = sigma2;
}

/* used for iClusterBayes, which only modle normal, poisson and binomial data */
void fillData3(dataType *dt,double *alpha,double *beta,double *con,int *cat,
	      double *sigma2,int *p,int *k,int *type){
  dt->type = type;
  dt->p = p;
  dt->k = k; 
  /* dt->c = c; */
  dt->alpha = alpha;
  dt->beta = beta;
  dt->con = con;
  dt->cat = cat;
  /* dt->class = class;
     dt->nclass = nclass; */
  dt->sigma2 = sigma2;
}

/*iCluster C programs -- implement Lasso, Enet, Group Lasso, Fussed Lasso, Group Fussed Lasso. 
  Last updated: 10/11/2010
  Qianxing Mo (qmo@bcm.edu), Dan L. Duncan Cancer Center, Baylor College of Medicine
*/


/*X is the t(X) of the original R code */
void  Mstep_glasso(int *p, int *k, int *n,double *B, double *X, double *Phivec, 
		   double *EZZt, double *EZ, double *lam, double *eps2){
  char *trans="N";
  int brow,bcol,xrow,xcol,ezrow, ezcol,ezztrow,j,incx,incy;
  double **bm, **xm, *IM, *ezzt,  *tempv1, *tempv2, w, ONE,ZERO; /* *phiv,  */
  int *IPIV, INFO;

  incx = 1;
  incy = 1;
  ONE = 1;
  ZERO = 0;

  brow = *p;
  bcol = *k;
  xrow = *p;
  xcol = *n;
  ezrow = *k;
  ezcol = *n;
  ezztrow = *k;

  IPIV = ivec(ezztrow);
  /*  IPIV = (int *)R_alloc(ezztrow, sizeof(int)); */
  bm = drowm(brow,bcol);
  xm = drowm(xrow,xcol);
  ezzt = dvec(ezztrow*ezztrow);
  /*  w = dvec(ezztrow); */

  IM = dvec(ezztrow*ezztrow);
  tempv1 = dvec(ezrow);
  tempv2 = dvec(ezztrow);

  dvtom(bm,B,brow,bcol);
  dvtom(xm,X,xrow,xcol);
  /* printvec(B,20); */
  /*  printmatrix(bm,10,2);   */
  /*  printmatrix(xm,10,8); */
  for(j=0; j< (brow); j++){
    /* j=0;
    dvcopy(w,lam,ezztrow);
    dvscale(w,ezztrow,Phivec[j]/sqrt(vecsum2(bm[j],bcol))); */
    w = Phivec[j]*(*lam)/sqrt(vecsum2(bm[j],bcol)); 
    /*   printf("-w-\n");
	 printvec(w,ezztrow); */

    dvcopy(ezzt,EZZt,ezztrow*ezztrow);
    diagplus(ezzt,ezztrow,w);
    diagm(IM,ezztrow,1);
    /*
    printf("-EZZt-\n");
    printvec(ezzt,ezztrow*ezztrow);
    */
    F77_CALL(dgesv)(&ezztrow,&ezztrow,ezzt,&ezztrow,IPIV,IM,&ezztrow,&INFO);
    F77_CALL(dgemv)(trans,&ezrow,&ezcol,&ONE,EZ,&ezrow, xm[j],&incx, &ZERO,tempv1, &incy FCONE);  
    F77_CALL(dgemv)(trans,&ezztrow,&ezztrow,&ONE,IM,&ezztrow,tempv1,&incx,&ZERO,tempv2,&incy FCONE); 
    /*   printvec(tempv2,ezztrow); */
    dvcopy(bm[j],tempv2,ezztrow);
  } 
  /*  printmatrix(bm,5,2); */
  editm(bm,(brow),(bcol),(*eps2));
  dmtov(B,bm,(brow),(bcol));
  /*  printf("- Good end -\n"); */
  /* printmatrix(bm,(brow),(bcol)); */
  /* dvcopy(B,bv,(bcol)*(brow)); */
  dmfree(bm,brow);
  dmfree(xm,xrow);
  Free(IM);
  Free(ezzt);
  /*  Free(phiv); */
  Free(tempv1);
  Free(tempv2);
  Free(IPIV);
}

/*X is the t(X) of the original R code */
void  Mstep_lasso(int *p, int *k, int *n,double *B, double *X, double *Phivec, 
		   double *EZZt, double *EZ, double *lam, double *eps2){
  char *trans = "N";
  int brow,bcol,xrow,xcol,ezrow, ezcol,ezztrow,j,incx,incy;
  double **bm, **xm, *IM, *ezzt, *tempv1, *tempv2, *w, ONE,ZERO;
  int *IPIV, INFO;

  incx = 1;
  incy = 1;
  ONE = 1;
  ZERO = 0;

  brow = *p;
  bcol = *k;
  xrow = *p;
  xcol = *n;
  ezrow = *k;
  ezcol = *n;
  ezztrow = *k;

  IPIV = ivec(ezztrow); 
  bm = drowm(brow,bcol);
  xm = drowm(xrow,xcol);
  ezzt = dvec(ezztrow*ezztrow);
  w = dvec(ezztrow);
 
  IM = dvec(ezztrow*ezztrow);
  tempv1 = dvec(ezrow);
  tempv2 = dvec(ezztrow);

  dvtom(bm,B,brow,bcol);
  dvtom(xm,X,xrow,xcol);
  /* printvec(B,20); */
  /*  printmatrix(bm,10,2);   */
  /*  printmatrix(xm,10,8); */
  for(j=0; j< (brow); j++){
    dvcopy(w,bm[j],ezztrow);
    fabsinv(w,w,ezztrow);
    dvscale(w,ezztrow,lam[0]*Phivec[j]);
    dvcopy(ezzt,EZZt,ezztrow*ezztrow);
    diagplusv(ezzt,ezztrow,w);
    /*  diagplus(ezzt,ezztrow,lam[1]); */
    diagm(IM,ezztrow,1);

    /*    printf("-EZZt-\n");
	  printvec(ezzt,ezztrow*ezztrow); */
    F77_CALL(dgesv)(&ezztrow,&ezztrow,ezzt,&ezztrow,IPIV,IM,&ezztrow,&INFO);
    F77_CALL(dgemv)(trans,&ezrow,&ezcol,&ONE,EZ,&ezrow, xm[j],&incx, &ZERO,tempv1, &incy FCONE);  
    F77_CALL(dgemv)(trans,&ezztrow,&ezztrow,&ONE,IM,&ezztrow,tempv1,&incx,&ZERO,tempv2,&incy FCONE); 
    /*  printvec(tempv2,ezztrow); */
    dvcopy(bm[j],tempv2,ezztrow);
  } 
/*  printf("eps2 = %f",*eps2);
    printmatrix(bm,5,2);  */
  editm(bm,(brow),(bcol),(*eps2));
  dmtov(B,bm,(brow),(bcol));
  /*  printmatrix(bm,(brow),(bcol)); */
  /*  dvcopy(B,bv,(bcol)*(brow));  */
  /*  printf("- good 6 -\n"); */
  dmfree(bm,brow);
  dmfree(xm,xrow);
  Free(IM);
  Free(ezzt);
  /*  Free(phiv); */
  Free(tempv1);
  Free(tempv2);
  Free(w);
  Free(IPIV);
}

/*X is the t(X) of the original R code */
void  Mstep_enet(int *p, int *k, int *n,double *B, double *X, double *Phivec, 
		   double *EZZt, double *EZ, double *lam, double *eps2){
  char *trans = "N";
  int brow,bcol,xrow,xcol,ezrow, ezcol,ezztrow,j,incx,incy;
  double **bm, **xm, *IM, *ezzt,  *tempv1, *tempv2, *w, ONE,ZERO; /* *phiv, */
  int *IPIV, INFO;

  incx = 1;
  incy = 1;
  ONE = 1;
  ZERO = 0;

  brow = *p;
  bcol = *k;
  xrow = *p;
  xcol = *n;
  ezrow = *k;
  ezcol = *n;
  ezztrow = *k;

  IPIV = ivec(ezztrow); 
  bm = drowm(brow,bcol);
  xm = drowm(xrow,xcol);
  ezzt = dvec(ezztrow*ezztrow);
  w = dvec(ezztrow);
 
  IM = dvec(ezztrow*ezztrow);
  tempv1 = dvec(ezrow);
  tempv2 = dvec(ezztrow);

  dvtom(bm,B,brow,bcol);
  dvtom(xm,X,xrow,xcol);
  /* printvec(B,20); */
  /*  printmatrix(bm,10,2);   */
  /*  printmatrix(xm,10,8); */
  for(j=0; j< (brow); j++){
    dvcopy(w,bm[j],ezztrow);
    fabsinv(w,w,ezztrow);
    dvscale(w,ezztrow,lam[0]*Phivec[j]);
    dvcopy(ezzt,EZZt,ezztrow*ezztrow);
    diagplusv(ezzt,ezztrow,w);
    diagplus(ezzt,ezztrow,lam[1]);
    diagm(IM,ezztrow,1);

    /*  printf("-EZZt-\n");
	printvec(ezzt,ezztrow*ezztrow); */
    F77_CALL(dgesv)(&ezztrow,&ezztrow,ezzt,&ezztrow,IPIV,IM,&ezztrow,&INFO);
    F77_CALL(dgemv)(trans,&ezrow,&ezcol,&ONE,EZ,&ezrow, xm[j],&incx, &ZERO,tempv1, &incy FCONE);  
    F77_CALL(dgemv)(trans,&ezztrow,&ezztrow,&ONE,IM,&ezztrow,tempv1,&incx,&ZERO,tempv2,&incy FCONE); 
    /*  printvec(tempv2,ezztrow); */
    dvcopy(bm[j],tempv2,ezztrow);
  } 
/*  printf("eps2 = %f",*eps2);
    printmatrix(bm,5,2); */
  editm(bm,(brow),(bcol),(*eps2));
  dmtov(B,bm,(brow),(bcol));
  /*  printmatrix(bm,(brow),(bcol)); */
  /*  dvcopy(B,bv,(bcol)*(brow));  */
  /*  printf("- good 6 -\n"); */
  dmfree(bm,brow);
  dmfree(xm,xrow);
  Free(IM);
  Free(ezzt);
  /*  Free(phiv); */
  Free(tempv1);
  Free(tempv2);
  Free(w);
  Free(IPIV);
}


void  Mstep_flasso(int *p, int *k,double *B, double *Phivec, double *EXZt,
		   double *EZZt, double *lam, double *eps2, int *id,int *idlen){
  char *trans = "N";
  int *IPIV;
  int brow,bcol,ezztrow,exztrow,exztcol,i,j,incx,incy,INFO,pk,pkpk,tempID1,tempID2,kj,colID,rowID;

  double **bm,*IM,*ezzt,**exzt,**M,**tempM, **qtilde,*W,*L,*MV,*tempm,*ctilde,*tempC,*bv,*tempv1,*tempv2;
  double ONE,ZERO;

  incx = 1;
  incy = 1;
  ONE = 1;
  ZERO = 0;

  brow = *p;
  bcol = *k;
  pk = brow*bcol;
  pkpk = pk*pk;
  ezztrow = *k;
  exztrow = *p;
  exztcol = *k;

  exzt = drowm(exztrow,exztcol);
  dvtom(exzt,EXZt,exztrow,exztcol);
  bm = drowm(brow,bcol);
  ezzt = dvec(ezztrow*ezztrow);
  M = drowm(pk,pk);

  IPIV = ivec(pk);
  IM = dvec(pkpk); /*identity matrix */
  tempv1 = dvec(bcol);
  bv = dvec(pk);
  W = dvec(pkpk);
  tempm = dvec(bcol*bcol); /* tempm = diag(1/tmp), tempm is a vector format of matrix; */
  dvtom(bm,B,brow,bcol);  /* bm is matrix B */
  dmtranv(bv,bm,brow,bcol); /* bv = t(bm) */

  for(i=0; i<pk; i++){
    if(bv[i]<0){
      bv[i] = -1.0*bv[i];
    }
  }

  dvinv(bv,pk,1.0);
  diagmv(W,pk,bv);

  /* printvec(B,20); */
  /* printmatrix(bm,10,2);   */
  /* printmatrix(xm,10,8); */

  for(j=0; j< (*idlen); j++){
    dvcopy(tempv1,bm[id[j]-1],bcol); /*note C is indexed from 0 */
    dvsub(tempv1,bm[id[j]],bcol);
    for(i=0; i<bcol; i++){
      tempv1[i] = fabs(tempv1[i]);
      if(tempv1[i] < (*eps2)){
	tempv1[i] = (*eps2);
      }
    }

    dvinv(tempv1,bcol,1.0);
    /* tempm = dvec(bcol*bcol); tempm = diag(1/tmp), tempm is a vector format of matrix; move this line to above in order to avoid warning */
    diagmv(tempm,bcol,tempv1);
    /* if(j==0){
           printvec(tempv1,bcol);
	   printvec(tempm,bcol*bcol); 
    } */
    tempID1 = 0; 
    kj = (*k)*id[j];
    for(colID = kj; colID <= (kj +(*k)-1); colID++){   
      for(rowID =(kj-(*k)); rowID <= (kj-1); rowID++){
 	M[rowID][colID] = tempm[tempID1];
	tempID1++;
      }
    }
  }

  tempM = drowm(pk,pk);
  dmtranm(tempM, M,pk,pk);
  dmadd(M,tempM,pk,pk);

  /* printmatrix(M,5,5); */
  tempv2 = dvec(pk); /*tempv2 = D  */
  dmrowsum(tempv2,M,pk,pk);

  L = dvec(pkpk);
  MV = dvec(pkpk);
  diagmv(L,pk,tempv2);
  dmtov(MV,M,pk,pk);
  dvsub(L,MV,pkpk);

  /*  printvec(L,10); */
  qtilde = drowm(pk,pk);
  ctilde = dvec(pk);
  tempC = dvec(exztcol);
  
  for(j=1; j<= (*p); j++){
    kj = (*k)*j;
    dvcopy(ezzt,EZZt,ezztrow*ezztrow);
    dvscale(ezzt,ezztrow,1.0/Phivec[j-1]); /* j is the R index */
    dvcopy(tempC,exzt[j-1],exztcol);
    dvscale(tempC,exztcol,1.0/Phivec[j-1]);
    tempID1 = 0;
    tempID2 = 0;
    for(colID = (kj-(*k)); colID <= (kj-1); colID++){  
      ctilde[colID] = tempC[tempID2];
      tempID2++;
      for(rowID =(kj-(*k)); rowID <= (kj-1); rowID++){
 	qtilde[rowID][colID] = ezzt[tempID1];
	tempID1++;
      }
    }  
  }
  /*  printmatrix(qtilde,5,5); */
  dvscale(W,pkpk,lam[0]);
  dvscale(L,pkpk,lam[1]);
  dvadd(W,L,pkpk);  /*W = P */

  dmtov(MV,qtilde,pk,pk); /*MV = qtilde */

  dvadd(W,MV,pkpk);

  diagm(IM,pk,1);   /* IM is identity vector */
  /*  printf("-EZZt-\n");
      printvec(EZZt,ezztrow*ezztrow); */
  F77_CALL(dgesv)(&pk,&pk,W,&pk,IPIV,IM,&pk,&INFO);
  F77_CALL(dgemv)(trans,&pk,&pk,&ONE,IM,&pk,ctilde,&incx,&ZERO,tempv2,&incy FCONE); 

  for(j=0; j<pk; j++){
    if(fabs(tempv2[j]) < (*eps2)){
      tempv2[j] = (*eps2);
    }
  }
 
  tempID1 = 0;
  for(i=0; i<brow; i++){
    for(j=0;j<bcol; j++){
      bm[i][j] = tempv2[tempID1];
      tempID1++;
    }
  }

  dmtov(B,bm,brow,bcol);

  dmfree(bm,brow);
  Free(IM);
  Free(ezzt);
  dmfree(exzt,exztrow);
  dmfree(M,pk);
  dmfree(tempM,pk);
  dmfree(qtilde,pk);
  Free(W);
  Free(L);
  Free(MV);
  Free(tempm);
  Free(ctilde);
  Free(tempC);
  Free(bv);
  /*  Free(phiv); */
  Free(tempv1);
  Free(tempv2);
  Free(IPIV);
}

/* w = dvec(n);  eigen values
   z= dvec(n*n); eigen vectors
*/
void eigen(double *A, int *row, double *w, double *z){
  char *jobz ="V", *range="A", *uplo="U";
  int n, lda,il, iu, m, ldz, *isuppz, lwork,*iwork,liwork,info,i,j,id;
  double vl,vu,abstol,*work;
  double *tempw, *tempz;
  tempw = dvec(*row);
  tempz = dvec((*row)*(*row));

  n = *row;
  lda= *row;
  vl = 0;
  vu = 1000;
  il = 1;
  iu = 1000;
  abstol = 0.0000001;
  m = 0;
  ldz = n;
  isuppz = ivec(2*n);
  work = dvec(100*n);
  lwork = 100*n;
  iwork = ivec(100*n);
  liwork = 100*n;
  info = 0;
  F77_CALL(dsyevr)(jobz, range,uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,&m,tempw,tempz,&ldz,isuppz,
		 work, &lwork,iwork, &liwork,&info FCONE FCONE FCONE);

  /*sort the eigen values and eigen vectors */
  id = 0;
  for(i=n; i>=1; i--){
    w[n-i] = tempw[i-1];
    for(j=(i*n-n);j<=(i*n-1);j++){
      z[id] = tempz[j];
      id++;
    }
  }
  Free(tempw);
  Free(tempz);
  Free(work);
  Free(isuppz);
  Free(iwork);
}


void lyap(double *B, double *P, double *Q, double *C, int *m, int *n){
  int i,j;
  char *transa="N", *transb="N";
  double alpha,beta;
  double *U,*lamvec,*V,*invV, *muvec, *ctilde, *btilde, *tempm, *tempv;
  alpha = 1;
  beta = 0;
  U = dvec((*m)*(*m));
  lamvec = dvec(*m);
  V = dvec((*n)*(*n));
  invV =  dvec((*n)*(*n));
  muvec = dvec(*n);
  ctilde = dvec((*m)*(*n));
  btilde = dvec((*m)*(*n));
  tempv = dvec((*m)*(*n));
  tempm = dvec((*m)*(*m));

  eigen(P,m,lamvec,U);
  eigen(Q,n,muvec,V);
  /* 
  printf("- V -\n");
  printvec(V,(*n)*(*n));
  printf("- Eigen Vector -\n");
  printvec(lamvec,10);
  printf("- Eigen Values -\n");
  printvec(muvec,*n); 
  */
  invsqm2(tempm, U, m); /*Note U is NOT change on exit */
  F77_CALL(dgemm)(transa,transb,m,n,n,&alpha,C,m,V,n,&beta,tempv,m FCONE FCONE);
  /* printvec(tempv,10); */
  F77_CALL(dgemm)(transa,transb,m,n,m,&alpha,tempm,m,tempv,m,&beta,ctilde,m FCONE FCONE);
  /*
  printf("-Ctilde 1 -\n");
  printvec(ctilde, 10);
  */

  for(i=0; i<(*n); i++){
    for(j=0; j<(*m); j++){
      btilde[(*m)*i+j] = ctilde[(*m)*i+j]/(lamvec[j]+muvec[i]);
    }
  }
  invsqm(invV, V, n);
  /*
  printf(" - btilde - \n");
  printvec(btilde,10);
  printf("- V -\n");
  printvec(invV,(*n)*(*n));
  */
  F77_CALL(dgemm)(transa,transb,m,n,n,&alpha,btilde,m,invV,n,&beta,tempv,m FCONE FCONE);
  F77_CALL(dgemm)(transa,transb,m,n,m,&alpha,U,m,tempv,m,&beta,B,m FCONE FCONE);
  Free(U);
  Free(lamvec);
  Free(V);
  Free(invV);
  Free(muvec);
  Free(ctilde);
  Free(btilde);
  Free(tempm);
  Free(tempv);
}


void  Mstep_gflasso(int *p,int *k,double *B,double *Phivec,double *EXZt,
		   double *EZZt,double *lam,double *eps2,int *id,int *idlen){

  char *transa = "N", *transb="N";;
  int brow,bcol,i,j,pp;
  double **bm, **M, *W,*L,*MV,*phi,*tempv1,*tempv2;
  double sumb, alpha, beta;

  alpha = 1;
  beta = 0;

  brow = *p;
  bcol = *k;
  pp = (*p)*(*p);

  bm = drowm(brow,bcol);
  M = drowm(*p,*p);
  W = dvec(pp);
  tempv1 = dvec(bcol);

  dvtom(bm,B,brow,bcol);  /* bm is matrix B */

  /*B^2 */
  tempv2 = dvec(*p);
  for(i=0; i<brow; i++){
    tempv2[i] = 0;
    for(j=0; j<bcol; j++){
      tempv2[i] = tempv2[i] + bm[i][j]*bm[i][j];
    }
    tempv2[i] = sqrt(tempv2[i]);
  }
  
  dvinv(tempv2, *p, 1.0);
  diagmv(W,brow,tempv2);

  for(j=0; j< (*idlen); j++){
    dvcopy(tempv1,bm[id[j]-1],bcol); /*note C is indexed from 0 */
    dvsub(tempv1,bm[id[j]],bcol);
    sumb = 0;
    for(i=0; i<bcol; i++){
      sumb = sumb + tempv1[i]*tempv1[i];
    }
    sumb = sqrt(sumb);
    if(sumb < (*eps2)){
      sumb = (*eps2);
    }
    M[id[j]-1][id[j]] = 1.0/sumb;
    M[id[j]][id[j]-1] =  M[id[j]-1][id[j]];
  }

  dmrowsum(tempv2,M,*p,*p);
 
  L = dvec(pp);
  MV = dvec(pp);
  diagmv(L,*p,tempv2);

  dmtov(MV,M,*p,*p);
 
  dvsub(L,MV,pp);
  /*
  printf("- L - \n");
  printvec(L,10);
  */
  dvscale(W,pp,lam[0]);
  dvscale(L,pp,lam[1]);
  dvadd(W,L,pp);  /*W = P */

  phi = dvec(pp);
  diagmv(phi,*p,Phivec);
  F77_CALL(dgemm)(transa,transb,p,p,p,&alpha,phi,p,W,p,&beta,MV,p FCONE FCONE);
  /*
  printf("- P -\n");
  printvec(MV,108);
  */
  lyap(B, MV,EZZt,EXZt,p,k);
  for(i=0; i<(brow*bcol); i++){
    if(fabs(B[i])<(*eps2)){
	B[i] = (*eps2);
    }
  } 

  dmfree(bm,brow);
  dmfree(M,*p);
  Free(W);
  Free(L);
  Free(MV);
  Free(phi);
  Free(tempv1);
  Free(tempv2);
  
}


/* dim(B) = p,k; dim(X) = n,p */
void iClusterCore(int *p, int *k, int *n, double *xtxdiag, double *X,double *B,double *EZ,double *EZZt,
		  double *Phivec, double *dif, int *iter, int *pvec,
		  double *lambda,double *eps,double *eps2, int *maxiter, int *lenT,int *method, int *ID, 
		  int *lenID){

  char *transN="N",*transT="T";
  int i, j,kk,pk,s,t; /* pp, */
  int ij;
  double *btp,*btpb, *EXZt,*tempX,*tempm0,*tempm1,*tempm2,*BOld,*PhivecOld, *XtXdiag; /* *tempm3, */
  double lbd1,lbd4,*lbd2,*lbd3,*lbd5;
  double alpha, beta,absdif;
  double ninv;

  double *xlist[*lenT],*phi[*lenT], *tempB[*lenT],*tempEXZt[*lenT],**xm, **xmt;

  /*  printf("lenT = %d %d",lenT[0],iter[0]); */
  alpha = 1;
  beta = 0;
  ninv = 1.0/(*n);

  kk = (*k)*(*k);
  pk = (*p)*(*k);
  /* kn = (*k)*(*n); */
  /* pp = (*p)*(*p); */

  btp = dvec(pk);
  btpb = dvec(kk);
  EXZt = dvec(pk);
  tempm0 = dvec(*p);
  tempm1 = dvec(kk);
  tempm2 = dvec(pk);
  // tempm3 = dvec(pp);
  BOld = dvec(pk);
  PhivecOld = dvec(*p);
  XtXdiag = dvec(*p);

  dvcopy(BOld,B,pk);
  dvcopy(PhivecOld,Phivec,*p);

  lbd2 = dvec(2);
  lbd3 = dvec(2);
  lbd5 = dvec(2);

  lbd1 = lambda[0];        /* lasso */
  lbd2[0]=lambda[1];       /* Enet  */
  lbd2[1] = lambda[2];
  lbd3[0]=lambda[3];      /* flasso */
  lbd3[1] = lambda[4];
  lbd4 = lambda[5];       /* glasso */
  lbd5[0] = lambda[6];   /* gflasso */
  lbd5[1] = lambda[7];

  s=0;
  for(t=0; t<(*lenT); t++){
    xlist[t] = dvec((*n)*pvec[t]);
    xm = drowm(*n,pvec[t]);
    xmt = drowm(pvec[t],*n);
    tempX = dvec((*n)*pvec[t]);
    dmsect(tempX,X,*n,0,(*n-1),s,(s+pvec[t]-1));
    dvtom(xm,tempX,*n,pvec[t]);
    dmtranm(xmt,xm,*n,pvec[t]);
    dmtov(xlist[t],xmt,pvec[t],*n);
    s = s+pvec[t];
    if(t==0){
      /*   printvec(xlist[t],pvec[t]*(*n)); */
    }
    dmfree(xm,*n);
    dmfree(xmt,pvec[t]);
    Free(tempX);
    phi[t] = dvec(pvec[t]);
    tempB[t] = dvec(pvec[t]*(*k));
    tempEXZt[t] = dvec(pvec[t]*(*k));
    /*    printf("- %d - \n",t); */
  }
  
  while((*dif) > (*eps) && (*iter) <(*maxiter)){ 
    *iter = *iter + 1;
    /* removed "dvcopy(btp,B,pk);" as it is redundant */
    /*   printf("- iter %d \n", *iter); */
    for(i=0; i< (*k); i++){
      for(j=0; j< (*p); j++){
	ij = i*(*p)+j; /* ith row, jth column index */
	/* btp stores Phi^{-1}B; so need to divide not multiply */
	btp[ij] = B[ij]/Phivec[j];
      }
    }
    /* btp stores Phi^{-1}B; transpose and multiply by B */
    F77_CALL(dgemm)(transT,transN,k,k,p,&alpha,btp,p,B,p,&beta,btpb,k FCONE FCONE);
 
    diagplus(btpb,*k,1);
    if((*iter)==0){
      /*     printf("- btpb -\n");
	     printvec(btpb,kk); */
    }
 
    invsqm(tempm1,btpb,k);  /* solve(BtinvPhiB); BtinvPhiB is not used anymore,so use invsqm  */
    if((*iter)==0){
      /*   printf("- inv btpb -\n");
	   printvec(tempm1,kk); */
    }
   
    F77_CALL(dgemm)(transN,transN,p,k,k,&alpha,btp,p,tempm1,k,&beta,tempm2,p FCONE FCONE); /*temp2 = tmp */
    if(*iter==0){
      /*   printf("- tmp -\n");
	   printvec(tempm2, 528); */
    }
   
    F77_CALL(dgemm)(transT,transT,k,n,p,&alpha,tempm2,p,X,n,&beta,EZ,k FCONE FCONE); /* EZ = t(tmp)%*%t(X) */
    if((*iter)==0){
      /*   printf("- EZ -\n");
	   printvec(EZ, 10); */
    }
    
    F77_CALL(dgemm)(transT,transT,p,k,n,&alpha,X,n,EZ,k,&beta,EXZt,p FCONE FCONE); /* EZZt = t(X)%*%t(EZ) */
    F77_CALL(dgemm)(transT,transN,k,k,p,&alpha,B,p,tempm2,p,&beta,EZZt,k FCONE FCONE);
    dvscale(EZZt,kk,-1.0);
    diagplus(EZZt,*k,1);
    dvscale(EZZt,kk,*n);
    F77_CALL(dgemm)(transN,transT,k,k,n,&alpha,EZ,k,EZ,k,&beta,tempm1,k FCONE FCONE);
    dvadd(EZZt,tempm1,kk);
    if((*iter)==0){
      /*   printf("- EZZt -\n");
	   printvec(EZZt, kk); */
    }   

    s=0;
    for(t=0; t<(*lenT); t++){
      /*    printf("s = %d, %d\n",s, s+pvec[t]-1); */
      dmsect(tempB[t],B,*p,s,(s+pvec[t]-1),0,(*k-1));
      dvsect(phi[t],Phivec,s,(s+pvec[t]-1));
      if(method[t] == 1){
	Mstep_lasso(&(pvec[t]),k,n,tempB[t],xlist[t],phi[t],EZZt,EZ,&lbd1,eps2);
	/*	if(*iter == 1 && t==0){
	  printvec(tempB[t],pvec[t]*(*k));
	  } */
	dmreplace(B,tempB[t],*p,s,(s+pvec[t]-1),0,(*k-1));
      }else if(method[t] == 2){
	Mstep_enet(&(pvec[t]),k,n,tempB[t],xlist[t],phi[t],EZZt,EZ,lbd2,eps2);
	dmreplace(B,tempB[t],*p,s,(s+pvec[t]-1),0,(*k-1)); 
      }else if(method[t] == 3){
	dmsect(tempEXZt[t],EXZt,*p,s,(s+pvec[t]-1),0,(*k-1));
	Mstep_flasso(&(pvec[t]),k,tempB[t],phi[t],tempEXZt[t],EZZt,lbd3,eps2,ID,lenID); 
	dmreplace(B,tempB[t],*p,s,(s+pvec[t]-1),0,(*k-1));
      }else if(method[t] == 4){
	Mstep_glasso(&(pvec[t]),k,n,tempB[t],xlist[t],phi[t],EZZt,EZ,&lbd4,eps2);
	dmreplace(B,tempB[t],*p,s,(s+pvec[t]-1),0,(*k-1));
      }else if(method[t] == 5){
	dmsect(tempEXZt[t],EXZt,*p,s,(s+pvec[t]-1),0,(*k-1));
	Mstep_gflasso(&(pvec[t]),k,tempB[t],phi[t],tempEXZt[t],EZZt,lbd5,eps2,ID,lenID); 
	dmreplace(B,tempB[t],*p,s,(s+pvec[t]-1),0,(*k-1));
      }
      /*  dmreplace(B,tempB[t],*p,s,(s+pvec[t]-1),0,(*k-1)); */
      s = s + pvec[t];
    }

    // this block is commented out since the computation can be done more efficiently  
    // F77_CALL(dgemm)(transN,transT,p,p,k,&alpha,EXZt,p,B,p,&beta,tempm3,p);
    // diagv(tempm0,tempm3,*p);
    // 
    // /* diagonal of t(X) %*% X */
    // dvcopy(XtXdiag,xtxdiag,*p);
    // /* diagonal of t(X) %*% X - EXZt%*%t(B) */
    // dvsub(XtXdiag,tempm0,*p);
    // /* diagonal of t(X) %*% X - EXZt%*%t(B) - B %*% t(EXZt) */
    // dvsub(XtXdiag,tempm0,*p);
    // 
    // F77_CALL(dgemm)(transN,transN,p,k,k,&alpha,B,p,EZZt,k,&beta,tempm2,p); /* B%*%EZZt */
    // F77_CALL(dgemm)(transN,transT,p,p,k,&alpha,tempm2,p,B,p,&beta,tempm3,p);
    // diagv(tempm0,tempm3,*p);
    // dvadd(XtXdiag,tempm0,*p);
    // dvscale(XtXdiag,*p,1.0/(*n));
    // dvcopy(Phivec,XtXdiag,*p);
    // /*    
    // if((*iter)==0){
    //   printf("- Phivec -\n");
    //   printvec(XtXdiag, 20);
    // }      
    
    /* direct computation of Phivec; don't need matrix product for just the diagonal */
    F77_CALL(dgemm)(transN,transN,p,k,k,&alpha,B,p,EZZt,k,&beta,tempm2,p FCONE FCONE); /* B%*%EZZt */
    for(j=0; j< (*p); j++){             /* j is the row index */
      Phivec[j] = xtxdiag[j];           /* initialize Phivec  */
      for(i=0; i< (*k); i++){           /* i the column index */
	ij = i*(*p)+j;                  /* jth row ith col index */
        /* inner product of jth row of (B%*%EZZt - 2*X%*%EZ) and jt row of B */
	Phivec[j] = Phivec[j] + (tempm2[ij] - 2.0*EXZt[ij])*B[ij];
      }
      /* divide by n */
      Phivec[j] = Phivec[j]*ninv;
    }

    *dif = 0;
    
    for(i=0; i<pk; i++){
      absdif = fabs(B[i]-BOld[i]);
      if(absdif > (*dif)){
	*dif = absdif;
      }
    }

    for(i=0; i<(*p); i++){
      absdif = fabs(Phivec[i]-PhivecOld[i]);
      if(absdif > (*dif)){
	*dif = absdif;
      }
    }

    dvcopy(BOld,B,pk);
    dvcopy(PhivecOld,Phivec,*p);
  }
  
  Free(btp);
  Free(btpb);
  Free(EXZt);
  Free(tempm0);
  Free(tempm1);
  Free(tempm2);
  // Free(tempm3);
  Free(BOld);
  Free(PhivecOld);
  Free(XtXdiag);
  Free(lbd2);
  Free(lbd3);
  Free(lbd5); 
  
  for(t=0; t<(*lenT); t++){
    Free(xlist[t]);
    Free(phi[t]);
    Free(tempB[t]);
    Free(tempEXZt[t]);
  }
  /*
  Free(xlist);
  Free(phi);
  Free(tempB);
  Free(tempEXZt);
  */
}

/*
SEXP invsqm(SEXP A, SEXP B){
  int i, n, *dimA, *IPIV, INFO;
  double *aptr, *bptr, *ansptr;
  SEXP ans;
  dimA = dim(A);
  n = dimA[0];
  PROTECT(ans = allocMatrix(REALSXP,dimA[0],dimA[1]));
  ansptr = REAL(ans);
  aptr = REAL(A);
  bptr = REAL(B);
  INFO = 0;
  IPIV = (int *)R_alloc(n, sizeof(int));
  F77_CALL(dgesv)(&n,&n,aptr,&n,IPIV,bptr,&n,&INFO);
  ansptr = aptr;
  return(ans);
}

SEXP testfun(SEXP A){
  int i,j,*dimA;
  double *pta, v=0;
  dimA = dim(A);
  PROTECT(A=coerceVector(A,REALSXP));
  pta = REAL(A);

  for(i=0; i< dimA[0]; i++){
    for(j=0; j < dimA[1]; j++){
      v += pta[i+j*dimA[0]];
      printf("%f ",pta[i + j*dimA[0]]);
    }
    printf("\n");
  }
  UNPROTECT(1);
  return R_NilValue;
}

*/


/* function for debuging */

/*
void lyap2(double *B, double *P, double *Q, double *C, double *UU, double *lam, double *VV, double *mu,
	   int *m, int *n){
  int i,j;
  char *transa="N", *transb="N";
  double alpha,beta;
  double *U,*lamvec,*V,*invV, *muvec, *ctilde, *btilde, *tempm, *tempv;
  alpha = 1;
  beta = 0;
  U = dvec((*m)*(*m));
  lamvec = dvec(*m);
  V = dvec((*n)*(*n));
  invV =  dvec((*n)*(*n));
  muvec = dvec(*n);
  ctilde = dvec((*m)*(*n));
  btilde = dvec((*m)*(*n));
  tempv = dvec((*m)*(*n));
  tempm = dvec((*m)*(*m));

  eigen(P,m,lamvec,U);
  eigen(Q,n,muvec,V);

  dvcopy(U,UU,(*m)*(*m));
  dvcopy(V,VV,(*n)*(*n));
  dvcopy(lamvec,lam,(*m));
  dvcopy(muvec,mu,(*n));
 
  printf("- U initial -\n");
  printvec(U,20);

  printf("- Eigen Vector -\n");
  printvec(lamvec,10);
  printf("- Eigen Values -\n");
  printvec(muvec,*n); 

  invsqm(tempm, U, m);
  printf("- U af invsqm -\n");
  printvec(U,20);
  F77_CALL(dgemm)(transa,transb,m,n,n,&alpha,C,m,V,n,&beta,tempv,m);
  printvec(tempv,10); 
  F77_CALL(dgemm)(transa,transb,m,n,m,&alpha,tempm,m,tempv,m,&beta,ctilde,m);
  
  printf("-Ctilde 1 -\n");
  printvec(ctilde, 10);
 
  for(i=0; i<(*n); i++){
    for(j=0; j<(*m); j++){
      btilde[(*m)*i+j] = ctilde[(*m)*i+j]/(lamvec[j]+muvec[i]);
    }
  }
  invsqm(invV, V, n);

  printf(" - btilde - \n");
  printvec(btilde,10);

  printf("- V -\n");
  printvec(invV,(*n)*(*n));

  F77_CALL(dgemm)(transa,transb,m,n,n,&alpha,btilde,m,invV,n,&beta,tempv,m);
  printf("- tempv -\n");
  printvec(tempv,(*m)*(*n)); 
  F77_CALL(dgemm)(transa,transb,m,n,m,&alpha,U,m,tempv,m,&beta,B,m);
  printf("- U -\n");
  printvec(U,20); 
  Free(U);
  Free(lamvec);
  Free(V);
  Free(invV);
  Free(muvec);
  Free(ctilde);
  Free(btilde);
  Free(tempm);
  Free(tempv);
}

*/
/*giCluster.c is the merged file of glmnetc.c and mixmcmc.c */
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

/* rsq = dev.ratio = (1-dev/nulldev) */
void elnetC(double *a0,double *beta,int *df,double *x,double *y,int *nobs,int *nvars,
	    double *alpha,double *lambda,double *rsq);

/* a0[p],beta[p][k],df[p],x[nobs][k],y[nobs][p], y = a0 + beta'x */
void elnetBatch(double *a0,double *beta,double *sigma2,double *x,double *y,int *nobs,int *k,int *p,
		double *alpha,double *lambda);

/* used to calculate sum of deviance and null deviance */
void elnetBatchDev(double *a0,double *beta,double *sigma2,double *x,double *y,int *nobs,int *k,int *p,
		double *alpha,double *lambda,double *sumdev0,double *sumdev);

/* dev0 is null deviance; dev = dev.ratio = (1-dev/nulldev) */
void fishnetC(double *a0,double *beta,int *df,double *x,double *y,int *nobs,int *nvars,
	      double *alpha,double *lambda,double *dev0,double *dev);

void fishnetBatch(double *a0,double *beta,double *x,double *y,int *nobs,int *k,int *p,
		  double *alpha,double *lambda);

/* used to calculate sum of deviance and null deviance */
void fishnetBatchDev(double *a0,double *beta,double *x,double *y,int *nobs,int *k,int *p,
		  double *alpha,double *lambda,double *sumdev,double *sumnulldev);

/* dev0 is null deviance; dev = dev.ratio = (1-dev/nulldev) */
void lognetC(double *a0,double *beta,int *df,double *x,int *y,int *nobs,int *nvars,
	     double *alpha,double *lambda,int *nclass,int *family,double *dev0,double *dev);

void lognetBatch(double *a0,double *beta,double *x,int *y,int *nobs,int *k,int *p,
		 double *alpha,double *lambda,int *nclass,int *maxclass,int *family);

/* used to calculate sum of deviance and null deviance */
void lognetBatchDev(double *a0,double *beta,double *x,int *y,int *nobs,int *k,int *p,
		 double *alpha,double *lambda,int *nclass,int *maxclass,int *family,
		 double *dev0,double *sumdev);

/* this function only work when lmu == 1 (nlam==1) */
/* get the coefficients from the elnet fortran program */
/* nin == maxnin when lambda is a scalar; number of compressed coefs for each solution */
/* beta is a vector with nvars elements; nx is not used here  */
void getbeta(double *beta,int *df,int *nin,int *nvars,int *ia,double *ca){
  int i;
  int *ja, *oja;
  double *tempca;

  ja = ivec(*nin);
  oja = ivec(*nin);
  tempca = dvec(*nin);

  *df = 0;
  for(i=0; i<(*nvars); i++){
    beta[i] = 0;
  }

  if((*nin)>0){
    for(i=0; i<(*nin); i++){
      ja[i] = ia[i];
      oja[i] = i;
      /*      tempca[i] = ca[i]; */
 
      if(ca[i] != 0){
	*df = *df + 1;
      }
    }
    
    R_qsort_int_I(ja,oja,1,*nin); /*ja and oja has been sorted in increasing order */

    for(i=0; i<(*nin); i++){
      tempca[i] = ca[oja[i]];
    }
    
    for(i=0; i<(*nin); i++){
      beta[ja[i]-1] = tempca[i]; /*because fortran is indexed from 1 */
    }
  }

  Free(ja);
  Free(oja);
  Free(tempca);
}

/* this function only work when lmu == 1 (nlam==1) */
/* get the coefficients from the lognet fortran program */
/* nin == maxnin when lambda is a scalar; number of compressed coefs for each solution */
/* beta is a matrix(nvars x nc); nx is not used here  */
void getbetaMult(double *beta,int *df,int *nin,int *nvars,int *nc,int *ia,double *ca){
  int i, j,id,nozero;
  int *ja, *oja;
  double **tempca1, **cacopy, **tempca2;

  ja = ivec(*nin);
  oja = ivec(*nin);
  tempca1 = drowm(*nin,*nc);
  cacopy = drowm(*nvars,*nc);
  tempca2 = drowm(*nvars,*nc);

  id = 0;
  for(j=0; j<(*nc); j++){
    for(i=0; i<(*nvars); i++){
      cacopy[i][j] = ca[id];
      id = id + 1;
    }
  }

  *df = 0;
  if((*nin)>0){
    for(i=0; i<(*nin); i++){
      ja[i] = ia[i];
      oja[i] = i;

      nozero = 0;
      for(j=0; j<(*nc); j++){
	if(cacopy[i][j] != 0){
	  nozero = 1;
	}
      }

      *df = *df + nozero;
    }
    
    R_qsort_int_I(ja,oja,1,*nin); /*ja and oja has been sorted in increasing order */

    for(i=0; i<(*nin); i++){
      for(j=0; j<(*nc); j++){
	tempca1[i][j] = cacopy[oja[i]][j];
      }
    }

    for(i=0; i<(*nin); i++){
      for(j=0; j<(*nc); j++){
	tempca2[ja[i]-1][j] = tempca1[i][j]; /*because fortran is indexed from 1 */
      }
    }
  }

  dmtov(beta,tempca2, *nvars, *nc);

  Free(ja);
  Free(oja);
  dmfree(tempca1, *nin);
  dmfree(cacopy, *nvars);
  dmfree(tempca2, *nvars);
}


/* elnetC call glmnet Fortran function elnet_; it work only when lambda is a scalar (nlam==1)  */
/* nobs = dim(x)[1], nvar=dim(x)[2], ulam=lambda */
/* vp = penalty.factor = rep(1,nvars),weights=rep(1,nobs) */
/* a0, beta[nvars],df,rsq are output; rsq is the dev.ratio = 1 - dev/nulldev */

void elnetC(double *a0,double *beta,int *df,double *x,double *y,int *nobs,int *nvars,
	    double *alpha,double *lambda,double *rsq){

  double *weights,*vp,*ca,alm=0,ulam,flmin=1.0,thresh=1e-7; /* meany; */
  int *ia,nin=0,i,ka,ne,nx,jd=0,nlam=1,isd=1,maxit=1000,lmu=0,nlp=0,jerr=0;

  *rsq = 0;
  *a0 = 0;
  ne = (*nvars) + 1; /* ne = dfmax + 1, maximum number of variables allowed to enter largest model */
  if(*nvars < 2*ne){  /* nx = dfmax = min(dfmax*2,nvars); nx will always == nvars if nvars > 0 */
    nx = *nvars;     /* maximum number of variables allowed to enter all models */
  }else{
    nx = 2*ne;
  }

  ca = dvec(nx*nlam);
  ia = ivec(nx);
  weights = dvec(*nobs);
  vp = dvec(*nvars);
  ulam = *lambda;

  for(i=0; i<(*nobs); i++){
    weights[i] = 1;
  }

  for(i=0; i<(*nvars); i++){
    vp[i] = 1; /* relative penalties for each predictor variable*/
  }

  if(*nvars < 500){
    ka = 1;  /* type.gaussian == covariance;  */
  }else{
    ka = 2;  /* type.gaussian == naive */
  }

  /*x, y, weights will be overwritten by elnet_ */
  elnet_(&ka,alpha,nobs,nvars,x,y,weights,&jd,vp,&ne,&nx,&nlam,&flmin,&ulam,&thresh,
	 &isd,&maxit,&lmu,a0,ca,ia,&nin,rsq,&alm,&nlp,&jerr);

  if(jerr == 0){
    getbeta(beta,df,&nin,nvars,ia,ca);
  }else if(jerr > 0){
    for(i=0; i<(*nvars); i++){
      beta[i] = 0;
    }
    *df = 0;
    Rprintf("Fatal Error! All beta values are set to zeros.");
  }else{
    for(i=0; i<(*nvars); i++){
      beta[i] = 0;
    }
    *df = 0;
    /*    Rprintf("Non Fatal Error! All beta values are set to zeros.");  */
  }

  Free(ca);
  Free(ia);
  Free(weights);
  Free(vp);
}


/* a0[p],beta[p][k],df[p],x[nobs][k],y[nobs][p], y = a0 + beta'x */
void elnetBatch(double *a0,double *beta,double *sigma2,double *x,double *y,int *nobs,int *k,int *p,
		double *alpha,double *lambda){

  int i, j, id, df,nk;
  double **bm, **xm, **ym,*xcopy, *ycopy,res,res2,rsq;

  nk = (*nobs)*(*k);
  bm = drowm(*p, *k);
  xm = drowm(*nobs, *k);
  ym = drowm(*p, *nobs); /* p x nobs */
  xcopy = dvec(nk);
  ycopy = dvec(*nobs);

  dvtom(xm, x, *nobs, *k);
  /* ym = t(y) */
  id = 0;
  for(i=0; i<(*p); i++){
    for(j=0; j<(*nobs); j++){
      ym[i][j] = y[id];
      id = id + 1;
    }
  }

  for(i=0; i<(*p); i++){
    /*x, y will be overwritten, thus have to copy  */
    dvcopy(xcopy,x,nk);
    dvcopy(ycopy,ym[i],*nobs);
    elnetC(&a0[i], bm[i],&df,xcopy,ycopy,nobs,k,alpha,lambda,&rsq); 

    res2 = 0;
    for(id=0; id<(*nobs); id++){
      res = ym[i][id] - a0[i] - dvdot(xm[id],bm[i],*k);
      res2 = res2 + res*res; 
    }
    sigma2[i] = res2/(*nobs-1-df);
  }

  dmtov(beta,bm,*p,*k);

  dmfree(bm, *p);
  dmfree(xm, *nobs);
  dmfree(ym, *p);
  Free(xcopy);
  Free(ycopy);
}

double nulldev(double *y,int n){
  int i;
  double meany, nd,dif;

  meany = 0;
  nd = 0;

  for(i=0; i<n; i++){
    meany = meany + y[i];
  }
  meany = meany/n;

  for(i=0; i<n; i++){
    dif = y[i] - meany;
    nd = nd + dif*dif;
  }
  return(nd);
} 


/* a0[p],beta[p][k],df[p],x[nobs][k],y[nobs][p], y = a0 + beta'x */
void elnetBatchDev(double *a0,double *beta,double *sigma2,double *x,double *y,int *nobs,int *k,int *p,
		   double *alpha,double *lambda,double *sumdev0,double *sumdev){

  int i, j, id, df,nk;
  double **bm, **xm, **ym,*xcopy, *ycopy,res,res2,rsq,dev0,dev;

  *sumdev0 = 0;
  *sumdev = 0;

  nk = (*nobs)*(*k);
  bm = drowm(*p, *k);
  xm = drowm(*nobs, *k);
  ym = drowm(*p, *nobs); /* p x nobs */
  xcopy = dvec(nk);
  ycopy = dvec(*nobs);

  dvtom(xm, x, *nobs, *k);
  /* ym = t(y) */
  id = 0;
  for(i=0; i<(*p); i++){
    for(j=0; j<(*nobs); j++){
      ym[i][j] = y[id];
      id = id + 1;
    }
  }

  for(i=0; i<(*p); i++){
    /*x, y will be overwritten, thus have to copy  */
    dvcopy(xcopy,x,nk);
    dvcopy(ycopy,ym[i],*nobs);
    elnetC(&a0[i], bm[i],&df,xcopy,ycopy,nobs,k,alpha,lambda,&rsq); 
  
    /* calculate deviance and null deviance */  
    dev0 = nulldev(ycopy,*nobs);
    dev = dev0*(1-rsq);
    *sumdev0 = *sumdev0 + dev0;
    *sumdev = *sumdev + dev;

    res2 = 0;
    for(id=0; id<(*nobs); id++){
      res = ym[i][id] - a0[i] - dvdot(xm[id],bm[i],*k);
      res2 = res2 + res*res; 
    }
    sigma2[i] = res2/(*nobs-1-df);
  }

  dmtov(beta,bm,*p,*k);

  dmfree(bm, *p);
  dmfree(xm, *nobs);
  dmfree(ym, *p);
  Free(xcopy);
  Free(ycopy);
}


/* fishnetC call glmnet Fortran function elnet_; it work only when lambda is a scalar (nlam==1)  */
/* nobs = dim(x)[1], nvar=dim(x)[2], ulam=lambda */
/* vp = penalty.factor = rep(1,nvars),weights=rep(1,nobs) */
/* a0, beta[nvars],df,dev0 (null deviance) and dev (dev.ratio) are output;  */

void fishnetC(double *a0,double *beta,int *df,double *x,double *y,int *nobs,int *nvars,
	      double *alpha,double *lambda,double *dev0,double *dev){
  double *weights,*offset,*vp,*ca,alm=0,ulam,flmin=1.0,thresh=1e-7;
  int *ia,i,ne,nx,nin=0,jd=0,nlam=1,isd=1,maxit=1000,lmu=0,nlp=0,jerr=0;

  *dev0 = 0;
  *dev = 0;
  ne = (*nvars) + 1; /* ne = dfmax + 1, maximum number of variables allowed to enter largest model */
  if(*nvars < 2*ne){  /* nx = dfmax = min(dfmax*2,nvars); nx will always == nvars if nvars > 0 */
    nx = *nvars;     /* maximum number of variables allowed to enter all models */
  }else{
    nx = 2*ne;
  }

  ca = dvec(nx*nlam);
  ia = ivec(nx);
  weights = dvec(*nobs);
  offset = dvec(*nobs);
  vp = dvec(*nvars);
  ulam = *lambda;

  for(i=0; i<(*nobs); i++){
    weights[i] = 1;
    offset[i] = 0;
  }

  for(i=0; i<(*nvars); i++){
    vp[i] = 1; /* relative penalties for each predictor variable*/
  }

  fishnet_(alpha,nobs,nvars,x,y,offset,weights,&jd,vp,&ne,&nx,&nlam,&flmin,&ulam,&thresh,
	   &isd,&maxit,&lmu,a0,ca,ia,&nin,dev0,dev,&alm,&nlp,&jerr);

  if(jerr == 0){
    getbeta(beta,df,&nin,nvars,ia,ca);
  }else if(jerr > 0){
    for(i=0; i<(*nvars); i++){
      beta[i] = 0;
    }
    *df = 0;
    Rprintf("Fatal Error! All beta values are set to zeros.");
  }else{
    for(i=0; i<(*nvars); i++){
      beta[i] = 0;
    }
    *df = 0;
    /*    Rprintf("Non Fatal Error! All beta values are set to zeros.");  */
  }

  Free(ca);
  Free(ia);
  Free(weights);
  Free(offset);
  Free(vp);
}

/* a0[p],beta[p][k],df[p],x[nobs][k],y[nobs][p], y = a0 + beta'x */
void fishnetBatch(double *a0,double *beta,double *x,double *y,int *nobs,int *k,int *p,
		  double *alpha,double *lambda){

  int i, j, id, df,nk;
  double **bm, **xm, **ym,*xcopy, *ycopy, dev0, dev;
  
  nk = (*nobs)*(*k);
  bm = drowm(*p, *k);
  xm = drowm(*nobs, *k);
  ym = drowm(*p, *nobs); /* p x nobs */
  xcopy = dvec(nk);
  ycopy = dvec(*nobs);

  dvtom(xm, x, *nobs, *k);

  id = 0;
  for(i=0; i<(*p); i++){
    for(j=0; j<(*nobs); j++){
      ym[i][j] = y[id];
      id = id + 1;
    }
  }

  for(i=0; i<(*p); i++){
    /*x, y will be overwritten, thus have to copy  */
    dvcopy(xcopy,x,nk);
    dvcopy(ycopy,ym[i],*nobs);
    fishnetC(&a0[i], bm[i],&df,xcopy,ycopy,nobs,k,alpha,lambda,&dev0,&dev);
  }

  dmtov(beta,bm,*p,*k);

  dmfree(bm, *p);
  dmfree(xm, *nobs);
  dmfree(ym, *p);
  Free(xcopy);
  Free(ycopy);
}

/* a0[p],beta[p][k],df[p],x[nobs][k],y[nobs][p], y = a0 + beta'x */
void fishnetBatchDev(double *a0,double *beta,double *x,double *y,int *nobs,int *k,int *p,
		     double *alpha,double *lambda,double *sumdev0,double *sumdev){

  int i, j, id, df,nk;
  double **bm, **xm, **ym,*xcopy, *ycopy, dev0, dev;
 
  *sumdev0 = 0;
  *sumdev = 0;

  nk = (*nobs)*(*k);
  bm = drowm(*p, *k);
  xm = drowm(*nobs, *k);
  ym = drowm(*p, *nobs); /* p x nobs */
  xcopy = dvec(nk);
  ycopy = dvec(*nobs);

  dvtom(xm, x, *nobs, *k);

  id = 0;
  for(i=0; i<(*p); i++){
    for(j=0; j<(*nobs); j++){
      ym[i][j] = y[id];
      id = id + 1;
    }
  }

  for(i=0; i<(*p); i++){
    /*x, y will be overwritten, thus have to copy  */
    dvcopy(xcopy,x,nk);
    dvcopy(ycopy,ym[i],*nobs);
    fishnetC(&a0[i], bm[i],&df,xcopy,ycopy,nobs,k,alpha,lambda,&dev0,&dev);
    *sumdev0 = (*sumdev0) + dev0;
    *sumdev = (*sumdev) + dev0*(1-dev);
  }

  dmtov(beta,bm,*p,*k);

  dmfree(bm, *p);
  dmfree(xm, *nobs);
  dmfree(ym, *p);
  Free(xcopy);
  Free(ycopy);
}


/* lognetC call glmnet Fortran function lognet_; it work only when lambda is a scalar (nlam==1)  */
/* y should be factor starting from 0; perform as.numeric(as.factor(y)) - 1 in R*/
/* nobs = dim(x)[1], nvar=dim(x)[2], ulam=lambda */
/* vp = penalty.factor = rep(1,nvars),weights=rep(1,nobs) */
/* a0, beta[nvars] and df are output;  */
/* if family = binomial; a0[1], beta[nvars]; if family = multinomial, a0[nclass],beta[nclass*nvars] */
/* nclass is a scalar */
void lognetC(double *a0,double *beta,int *df,double *x,int *y,int *nobs,int *nvars,
	     double *alpha,double *lambda,int *nclass,int *family,double *dev0,double *dev){
  double *ymat,*offset,*vp,*ca,alm=0,ulam,flmin=1.0,thresh=1e-7,a0mean=0;
  int *ia,nin=0,kopt,i,ne,nx,nc,jd=0,nlam=1,isd=1,maxit=1000,lmu=0,nlp=0,jerr=0;
  
  *dev0 = 0;
  *dev = 0;    
  if((*family) == 0){ /* binomial */
    nc = 1;
  }else{             /* if family !=0,even nclass == 2, treat y as multinomial case */
    nc = *nclass;
  }

  ne = (*nvars) + 1; /* ne = dfmax + 1, maximum number of variables allowed to enter largest model */
  if(*nvars < 2*ne){  /* nx = dfmax = min(dfmax*2,nvars); nx will always == nvars if nvars > 0 */
    nx = *nvars;     /* maximum number of variables allowed to enter all models */
  }else{
    nx = 2*ne;
  }

  kopt = 0;  /* use the exact Hessian */
  ulam = *lambda;

  ymat = dvec((*nobs)*(*nclass)); 
  ca = dvec(nx*nlam*nc);
  ia = ivec(nx);
  /*  weights = dvec(*nobs); */
  offset = dvec(*nobs*(*nclass));
  vp = dvec(*nvars);

  /* create a y matrix as lognet makes; y[i] must be the column index of ymat, 0:(nclass-1) */
  for(i=0; i< (*nobs); i++){
    ymat[i + y[i]*(*nobs)] = 1;
  }

  /*
  for(i=0; i<(*nobs); i++){
    weights[i] = 1;
    offset[i] = 0; 
  }
 */

  for(i=0; i<(*nvars); i++){
    vp[i] = 1; /* relative penalties for each predictor variable*/
  }

  /*  printf("alpha =%f, nobs=%d, nvars=%d, nc=%d,jd=%d,ne=%d,nx=%d,nlam=%d,flmin=%f,ulam=%f,thre=%f,isd=%d,maxit=%d,kopt=%d,lmu=%d\n", *alpha,*nobs,*nvars,nc,jd,ne,nx,nlam,flmin,ulam,thresh,isd,maxit,kopt,lmu); */
  /*  Rprintf("B/"); */
  lognet_(alpha,nobs,nvars,&nc,x,ymat,offset,&jd,vp,&ne,&nx,&nlam,&flmin,&ulam,&thresh,
	  &isd,&maxit,&kopt,&lmu,a0,ca,ia,&nin,dev0,dev,&alm,&nlp,&jerr);
  /*  Rprintf("A "); */
  if((*family)==0){
    if(jerr == 0){
      getbeta(beta,df,&nin,nvars,ia,ca);
    }else if(jerr > 0){
      for(i=0; i<(*nvars); i++){
	beta[i] = 0;
      }
      *df = 0;
      Rprintf("WARNING: Fatal Error! All beta values are set to zeros.\n");
    }else{
      for(i=0; i<(*nvars); i++){
	beta[i] = 0;
      }
      *df = 0;
      /*      Rprintf("Non Fatal Error! All beta values are set to zeros.\n");  */
    }
    /* reverse the sign for binomial data */    
    *a0 = *a0 *(-1);
    for(i=0; i<(*nvars); i++){
      if(beta[i] != 0){
	beta[i] = -beta[i];
      }
    }
  }else{
    a0mean = 0;
    for(i=0; i<(*nclass); i++){
      a0mean = a0mean + a0[i];
    }
    a0mean = a0mean/(*nclass);

    /* center a0 */
    for(i=0; i<(*nclass); i++){
      a0[i] = a0[i] - a0mean;
    }

    if(jerr == 0){
      getbetaMult(beta,df,&nin,nvars,nclass,ia,ca);
    }else if(jerr > 0){
      for(i=0; i<((*nvars)*(*nclass)); i++){
	beta[i] = 0;
      }
      *df = 0;
      Rprintf("Warning: Fatal Error! All beta values are set to zeros.\n");
    }else{
      for(i=0; i<((*nvars)*(*nclass)); i++){
	beta[i] = 0;
      }
      *df = 0;
      /*      Rprintf("Non Fatal Error! All beta values are set to zeros.\n");  */
    }
  }
    
  Free(ia);
  Free(ymat); 
  Free(ca);
  Free(offset);
  Free(vp);
}


/* a0[p],beta[p][k],df[p],x[nobs][k],y[nobs][p], y = a0 + beta'x */
/* nclass is vector, indicating the number of class for the columns of y */
/* maxclass is the maximum of nclass */
void lognetBatch(double *a0,double *beta,double *x,int *y,int *nobs,int *k,int *p,
		 double *alpha,double *lambda,int *nclass,int *maxclass,int *family){

  int i, j, id, df,nk, nc,nck, **ym, *ycopy;
  double **bm, **xm, *xcopy, **tempA0,dev0,dev;
  
  if((*family) == 0){ /* binomial */
    nc = 1;
  }else{
    nc = *maxclass;  /*assign the memory based on maxclass */
  }
  tempA0 = drowm(*p,nc);
  nck = (*k)*nc;
  nk = (*nobs)*(*k);
  bm = drowm(*p, nck);
  xm = drowm(*nobs, *k);
  ym = irowm(*p, *nobs); /* p x nobs */
  xcopy = dvec(nk);
  ycopy = ivec(*nobs);

  /*  dvtom(bm, beta, *p, *k); */
  dvtom(xm, x, *nobs, *k);

  id = 0;
  for(i=0; i<(*p); i++){
    for(j=0; j<(*nobs); j++){
      ym[i][j] = y[id];
      id = id + 1;
    }
  }
  /*  Rprintf(" --- lognetBatch before for loop .\n"); */
  for(i=0; i<(*p); i++){
    /*x, y will be overwritten, thus have to copy  */
    dvcopy(xcopy,x,nk);
    ivcopy(ycopy,ym[i],*nobs);
    /*    Rprintf("%d ",i); */
    lognetC(tempA0[i],bm[i],&df,xcopy,ycopy,nobs,k,alpha,lambda,&nclass[i],family,&dev0,&dev);
  }
  /*  Rprintf("\n --- lognetBatch after for loop .\n"); */
  dmtov(beta,bm,*p,nck);
  dmtov(a0,tempA0,*p,nc);

  dmfree(bm, *p);
  dmfree(xm, *nobs);
  imfree(ym, *p);
  dmfree(tempA0, *p);
  Free(xcopy);
  Free(ycopy);
}



/* a0[p],beta[p][k],df[p],x[nobs][k],y[nobs][p], y = a0 + beta'x */
/* nclass is vector, indicating the number of class for the columns of y */
/* maxclass is the maximum of nclass */
void lognetBatchDev(double *a0,double *beta,double *x,int *y,int *nobs,int *k,int *p,
		    double *alpha,double *lambda,int *nclass,int *maxclass,int *family,
		    double *sumdev0,double *sumdev){

  int i, j, id, df,nk, nc,nck, **ym, *ycopy;
  double **bm, **xm, *xcopy, **tempA0,dev0,dev;
  
  *sumdev0 = 0;
  *sumdev = 0;

  if((*family) == 0){ /* binomial */
    nc = 1;
  }else{
    nc = *maxclass;  /*assign the memory based on maxclass */
  }
  tempA0 = drowm(*p,nc);
  nck = (*k)*nc;
  nk = (*nobs)*(*k);
  bm = drowm(*p, nck);
  xm = drowm(*nobs, *k);
  ym = irowm(*p, *nobs); /* p x nobs */
  xcopy = dvec(nk);
  ycopy = ivec(*nobs);

  /*  dvtom(bm, beta, *p, *k); */
  dvtom(xm, x, *nobs, *k);

  id = 0;
  for(i=0; i<(*p); i++){
    for(j=0; j<(*nobs); j++){
      ym[i][j] = y[id];
      id = id + 1;
    }
  }

  for(i=0; i<(*p); i++){
    /*x, y will be overwritten, thus have to copy  */
    dvcopy(xcopy,x,nk);
    ivcopy(ycopy,ym[i],*nobs);
    lognetC(tempA0[i],bm[i],&df,xcopy,ycopy,nobs,k,alpha,lambda,&nclass[i],family,&dev0,&dev);
    *sumdev0 = (*sumdev0) + dev0;
    *sumdev = (*sumdev) + dev0*(1-dev);
  }

  dmtov(beta,bm,*p,nck);
  dmtov(a0,tempA0,*p,nc);

  dmfree(bm, *p);
  dmfree(xm, *nobs);
  imfree(ym, *p);
  dmfree(tempA0, *p);
  Free(xcopy);
  Free(ycopy);
}


/* the following code is from mixmcmc.c */
/* continous work on mixpca.c; handle mixed type of categorical data by modifying logMult,metroMix and mcmcXmix */ 
/* Qianxing Mo, Dan L. Duncan Cancer Center */
/*Last update: 7/5/2011  */


/*  logpnull <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
    logqnull <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta))) */
double logpnull(double etai){
  double res;
  if(etai < 0){
    res = etai -log1p(exp(etai));
  }else{
    res =  - log1p(exp(-etai));
  }
  return(res);
}

double logqnull(double etai){
  double res;
  if(etai < 0){
    res = -log1p(exp(etai));
  }else{
    res = -etai - log1p(exp(-etai));
  }
  return(res);
}

/* log Binomial likelihood for subject i */
/* alpha is p */
/* beta is p by k */
/* x is 1 by p  */
/* zi is 1 by k  */

double logBinom(double *zi,double *alpha,double *beta,int *xi,int *p,int *k){
  char *trans="N";
  int incx,incy,i;
  double ONE,loglike;
  double *eta;

  eta = dvec(*p);
  incx = 1;
  incy = 1;
  ONE = 1.0;
  
  dvcopy(eta,alpha,*p);
  F77_CALL(dgemv)(trans,p,k,&ONE,beta,p,zi,&incx,&ONE,eta,&incy FCONE);
  
  loglike = 0;
  for(i=0; i<(*p);i++){
    if(xi[i]==1){
      loglike = loglike + logpnull(eta[i]);
    }else{
      loglike = loglike + logqnull(eta[i]);
    }
  }
  
  /* The following line is for individual data set, not for mixed data set */
  /* de-comment this for testing single mcmcMult */ 
  /*  loglike = loglike - 0.5*F77_CALL(ddot)(k,zi,&incx,zi,&incy); */
  Free(eta);
  return(loglike);
}


/* log Binomial likelihood for subject i to n */
/* alpha is p */
/* beta is p by k */
/* X is n by p  */
/* Z is n by k  */
void logBinomAll(double *LogLike,double *Z,double *alpha,double *beta,int *X,int *n, int *p,int *k){
  char *trans="N";
  int incx,incy,i, j;
  double ONE,loglike;
  double *eta, **zm;
  int **xm;

  xm = irowm(*n,*p);
  ivtom(xm,X,(*n),(*p)); /* copy X to xm */
  zm = drowm(*n,*k);
  dvtom(zm,Z,*n,*k);    /* copy Z to zm */

  eta = dvec(*p);
  incx = 1;
  incy = 1;
  ONE = 1.0;
  
  loglike = 0;

  /* j is the index for subject */
  for(j=0; j<(*n); j++){
    dvcopy(eta,alpha,*p);
    F77_CALL(dgemv)(trans,p,k,&ONE,beta,p,zm[j],&incx,&ONE,eta,&incy FCONE);
    
    for(i=0; i<(*p);i++){
      if(xm[j][i]==1){
	loglike = loglike + logpnull(eta[i]);
      }else{
	loglike = loglike + logqnull(eta[i]);
      }
    }
  }

  *LogLike = loglike;

  /* The following line is for individual data set, not for mixed data set */
  /* de-comment this for testing single mcmcMult */ 
  /*  loglike = loglike - 0.5*F77_CALL(ddot)(k,zi,&incx,zi,&incy); */
  Free(eta);
  imfree(xm, *n);
  dmfree(zm, *n);
}


/* log Poisson likelihood for subject i */
/* alpha is p */
/* beta is p by k */
/* x is 1 by p  */
/* zi is 1 by k  */

double logPoisson(double *zi,double *alpha,double *beta,int *xi,int *p,int *k){
  char *trans="N";
  int incx,incy,i;
  double ONE,loglike;
  double *eta;

  eta = dvec(*p);
  incx = 1;
  incy = 1;
  ONE = 1.0;
  
  dvcopy(eta,alpha,*p);
  F77_CALL(dgemv)(trans,p,k,&ONE,beta,p,zi,&incx,&ONE,eta,&incy FCONE);
  
  loglike = 0;
  for(i=0; i<(*p);i++){
    loglike = loglike + xi[i]*eta[i] - exp(eta[i]);
  }

  /* The following line is for individual data set, not for mixed data set */
  /* de-comment this for testing single mcmcMult */  
  /*  loglike = loglike - 0.5*F77_CALL(ddot)(k,zi,&incx,zi,&incy); */
  /*  printf("%f \n",loglike); */
  Free(eta);
  return(loglike);
}

/* log Poisson likelihood for subject i to n */
/* alpha is p */
/* beta is p by k */
/* X is n by p  */
/* Z is n by k  */

void logPoissonAll(double *LogLike,double *Z,double *alpha,double *beta,int *X,int *n,int *p,int *k){
  char *trans="N";
  int incx,incy,i, j;
  double ONE,loglike;
  double *eta, **zm;
  int **xm;

  xm = irowm(*n,*p);
  ivtom(xm,X,(*n),(*p)); /* copy X to xm */
  zm = drowm(*n,*k);
  dvtom(zm,Z,*n,*k);    /* copy Z to zm */

  eta = dvec(*p);
  incx = 1;
  incy = 1;
  ONE = 1.0;

  loglike = 0;  
  /*loop over subject j, j = 1, ..., n */
  for(j=0; j<(*n); j++){
    dvcopy(eta,alpha,*p);
    F77_CALL(dgemv)(trans,p,k,&ONE,beta,p,zm[j],&incx,&ONE,eta,&incy FCONE);
  
    for(i=0; i<(*p);i++){
      loglike = loglike + xm[j][i]*eta[i] - exp(eta[i]);
    }
  }
  *LogLike = loglike;
  /* The following line is for individual data set, not for mixed data set */
  /* de-comment this for testing single mcmcMult */  
  /*  loglike = loglike - 0.5*F77_CALL(ddot)(k,zi,&incx,zi,&incy); */
  /*  printf("%f \n",loglike); */
  Free(eta);
  imfree(xm, *n);
  dmfree(zm, *n);
}

/* log normal likelihood for subject i */
/* alpha is p */
/* beta is p by k */
/* yi is 1 by p  */
/* zi is 1 by k  */
/* sigma2 is p */

double logNorm(double *zi,double *alpha,double *beta,double *sigma2,double *yi,int *p,int *k){
  char *trans="N";
  int incx,incy,i;
  double ONE,loglike;
  double *eta, *dif;

  eta = dvec(*p);
  dif = dvec(*p);
  incx = 1;
  incy = 1;
  ONE = 1.0;
  loglike = 0;

  /*Note, without copying alpha to eta because alpha will be substracted at the following code */
  F77_CALL(dgemv)(trans,p,k,&ONE,beta,p,zi,&incx,&ONE,eta,&incy FCONE);
  
  for(i=0;i<(*p); i++){
    dif[i] = yi[i] - alpha[i] - eta[i];
    loglike = loglike - log(sqrt(sigma2[i])) - 0.5*dif[i]*dif[i]/sigma2[i];
  }

  /* The following line is for individual data set, not for mixed data set */
  /* de-comment this for testing single mcmcMult */
  /*  loglike = loglike - 0.5*F77_CALL(ddot)(k,zi,&incx,zi,&incy); */
  /*  printf("%f \n",loglike); */
  Free(eta);
  Free(dif);
  return(loglike);
}


/* log normal likelihood for subject i to n */
/* alpha is p */
/* beta is p by k */
/* Y is n by p  */
/* Z is n by k  */
/* sigma2 is p */

void logNormAll(double *LogLike,double *Z,double *alpha,double *beta,double *sigma2,double *Y,int *n,int *p,int *k){
  char *trans="N";
  int incx,incy,i, j;
  double ONE,loglike;
  double *eta, *dif;
  double **ym, **zm;

  ym = drowm(*n,*p);
  dvtom(ym,Y,(*n),(*p)); /* copy X to xm */
  zm = drowm(*n,*k);
  dvtom(zm,Z,*n,*k);    /* copy Z to zm */

  eta = dvec(*p);
  dif = dvec(*p);
  incx = 1;
  incy = 1;
  ONE = 1.0;
  loglike = 0;

  /* loop over subject j, j = 1,..., n */
  for(j=0; j<(*n); j++){
    /*Note, without copying alpha to eta because alpha will be substracted at the following code */
    /* eta is a vector of zeros */
    F77_CALL(dgemv)(trans,p,k,&ONE,beta,p,zm[j],&incx,&ONE,eta,&incy FCONE);
  
    for(i=0;i<(*p); i++){
      dif[i] = ym[j][i] - alpha[i] - eta[i];
      loglike = loglike - log(sqrt(sigma2[i])) - 0.5*dif[i]*dif[i]/sigma2[i];
      eta[i] = 0; /* set eta to zero for F77_CALL(dgemv) */
    }
  }

  *LogLike = loglike;
  /* The following line is for individual data set, not for mixed data set */
  /* de-comment this for testing single mcmcMult */
  /*  loglike = loglike - 0.5*F77_CALL(ddot)(k,zi,&incx,zi,&incy); */
  /*  printf("%f \n",loglike); */
  Free(eta);
  Free(dif);
  dmfree(ym, *n);
  dmfree(zm, *n);
}


int indicator(int x, int y){
  int a;
  a = 0;
  if(x==y){
    a=1;
  }
  return(a);
}


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

double logMult(double *zi,double *alpha,double *beta,int *xi,int *class,int *nclass,int *p,int *C,int *K){
  /* char *trans="N"; */
  int incx,incy,j,c,k,betaID;
  double loglike, maxEtaj=0,expEtajc;
  double *etaj;
  double ** mAlpha; 
  double *betajc;

  betajc = dvec(*K);
  mAlpha = drowm(*p, *C);
  etaj = dvec(*C);
  incx = 1;
  incy = 1;
  /*  ONE = 1.0; */

  dvtom(mAlpha,alpha,*p,*C); 
  
  loglike = 0;
  betaID = 0;
  for(j=0; j<(*p);j++){
    for(c=0; c<(*C); c++){ /*important to loop over C because beta is p by k*C matrix */
      for(k=0; k<(*K); k++){
	betajc[k] = beta[betaID];
	betaID = betaID + 1;
      }
      if(c < nclass[j]){ /* nclass[j] may < *C; only calculate the necessary class */
	etaj[c] = mAlpha[j][c] + F77_CALL(ddot)(K,zi,&incx,betajc,&incy);
	if(c>0){
	  if(etaj[c] > maxEtaj){
	    maxEtaj = etaj[c];
	  }
	}else{
	  maxEtaj = etaj[c];
	} 
      }
    }
    
    expEtajc = 0;
    for(c=0; c<(nclass[j]); c++){ /* Only loop over nclass[j] because nclass may not equal to *C */
      expEtajc  += exp(etaj[c] - maxEtaj);
    }
    /*   expEtajc = log(expEtajc); */
    
    /* Only loop over nclass[j] because nclass[j] may < *C */
    for(c=0; c<(nclass[j]); c++){
      /*      loglike += indicator(xi[j],class[c]) * (etaj[c] -maxEtaj - log(expEtajc)); */
      if(xi[j] == class[c]){
	loglike += etaj[c] - maxEtaj - log(expEtajc);
      }
    }
  }

  /* The following line is for individual data set, not for mixed data set */
  /* de-comment this for testing single mcmcMult */
 
  /* loglike = loglike - 0.5*F77_CALL(ddot)(K,zi,&incx,zi,&incy); */
  /*  printf("%f \n", loglike); */
  dmfree(mAlpha,*p);
  Free(etaj);
  Free(betajc);
  return(loglike);
}



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
void logMultAll(double *LogLike,double *Z,double *alpha,double *beta,int *X,int *class,int *nclass,int *n,int *p,int *C,int *K){
  /* char *trans="N"; */
  int incx,incy,j,c,k,betaID,s;
  double loglike, maxEtaj=0,expEtajc;
  double *etaj, *betajc;;
  double **mAlpha, **zm; 
  int **xm;

  xm = irowm(*n,*p);
  ivtom(xm,X,(*n),(*p)); /* copy X to xm */
  zm = drowm(*n,*K);
  dvtom(zm,Z,*n,*K);    /* copy Z to zm */

  betajc = dvec(*K);
  mAlpha = drowm(*p, *C);
  etaj = dvec(*C);
  incx = 1;
  incy = 1;
  /*  ONE = 1.0; */
 
  dvtom(mAlpha,alpha,*p,*C); /* Convert alpha to matrix */
  
  loglike = 0;

  /* Loop over subject 1 to n */
  for(s=0; s<(*n); s++){
    betaID = 0;
    for(j=0; j<(*p);j++){
      for(c=0; c<(*C); c++){ /*important to loop over C because beta is p by k*C matrix */
	for(k=0; k<(*K); k++){
	  betajc[k] = beta[betaID];
	  betaID = betaID + 1;
	}
	if(c < nclass[j]){ /* nclass[j] may < *C; only calculate the necessary class */
	  etaj[c] = mAlpha[j][c] + F77_CALL(ddot)(K,zm[s],&incx,betajc,&incy);
	  if(c>0){
	    if(etaj[c] > maxEtaj){
	      maxEtaj = etaj[c];
	    }
	  }else{
	    maxEtaj = etaj[c];
	  } 
	}
      }

      expEtajc = 0;
      for(c=0; c<(nclass[j]); c++){ /* Only loop over nclass[j] because nclass may not equal to *C */
	expEtajc  += exp(etaj[c] - maxEtaj);
      }
      
      /* Only loop over nclass[j] because nclass[j] may < *C */
      for(c=0; c<(nclass[j]); c++){
	if(xm[s][j] == class[c]){
	  loglike += etaj[c] - maxEtaj - log(expEtajc);
	}
      }
    }
  }

  *LogLike = loglike;
  /* The following line is for individual data set, not for mixed data set */
  /* de-comment this for testing single mcmcMult */
 
  /* loglike = loglike - 0.5*F77_CALL(ddot)(K,zi,&incx,zi,&incy); */
  /*  printf("%f \n", loglike); */
  dmfree(mAlpha,*p);
  imfree(xm, *n);
  dmfree(zm, *n);
  Free(etaj);
  Free(betajc);
}

void metroBinom(double *zi,double *alpha,double *beta,int *xi,int *p,int *k,double *sigma,double *newz){
  double dif, *tryz;
  int i;
  tryz = dvec(*k);
  for(i=0; i<(*k); i++){
    tryz[i] = zi[i] + rnorm(0, *sigma);
  }
  dif = logBinom(tryz,alpha,beta,xi,p,k) - logBinom(zi,alpha,beta,xi,p,k);
  if((dif > 0) | (runif(0,1) < exp(dif))){
    dvcopy(newz, tryz, *k);
  }else{
    dvcopy(newz, zi, *k);
  }
  Free(tryz);
}


void metroPoisson(double *zi,double *alpha,double *beta,int *xi,int *p,int *k,double *sigma,double *newz){
  double dif, *tryz;
  int i;
  tryz = dvec(*k);
  for(i=0; i<(*k); i++){
    tryz[i] = zi[i] + rnorm(0, *sigma);
  }
  dif = logPoisson(tryz,alpha,beta,xi,p,k) - logPoisson(zi,alpha,beta,xi,p,k);
  if((dif > 0) | (runif(0,1) < exp(dif))){
    dvcopy(newz, tryz, *k);
  }else{
    dvcopy(newz, zi, *k);
  }
  Free(tryz);
}

/* beta = p  x (k x C) is the format in R; beta must be transposed using t(beta) before 
   passing to this function */
void metroMult(double *zi,double *alpha,double *beta,int *xi,int *class,int *nclass,int *p,int *c,int *k,double *sigma,double *newz){
  double dif, *tryz;
  int i;
  tryz = dvec(*k);
  for(i=0; i<(*k); i++){
    tryz[i] = zi[i] + rnorm(0, *sigma);
  }
  dif = logMult(tryz,alpha,beta,xi,class,nclass,p,c,k) - logMult(zi,alpha,beta,xi,class,nclass,p,c,k);
  if((dif > 0) | (runif(0,1) < exp(dif))){
    dvcopy(newz, tryz, *k);
  }else{
    dvcopy(newz, zi, *k);
  }
  Free(tryz);
}

void metroNormal(double *zi,double *mu,double *beta, double *sigma2,double *yi,int *p,int *k,double *sigma,double *newz){
  double dif, *tryz;
  int i;

  tryz = dvec(*k);
  for(i=0; i<(*k); i++){
    tryz[i] = zi[i] + rnorm(0, *sigma);
  }
  dif = logNorm(tryz,mu,beta,sigma2,yi,p,k) - logNorm(zi,mu,beta,sigma2,yi,p,k);
  if((dif > 0) | (runif(0,1) < exp(dif))){
    dvcopy(newz, tryz, *k);
  }else{
    dvcopy(newz, zi, *k);
  }
  Free(tryz);
}

/*if tryz is accepted, accept = accept + 1 */
void metroMix(double *zi,dataType *dt1,dataType *dt2,dataType *dt3,dataType *dt4, double *sigma,
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

  dif = 0;
  if((*ndt) > 0){
    switch(*(dt1->type)){
    case 1 : dif += logNorm(tryz,dt1->alpha,dt1->beta,dt1->sigma2,dt1->con,dt1->p,dt1->k) 
	- logNorm(zi,dt1->alpha,dt1->beta,dt1->sigma2,dt1->con,dt1->p,dt1->k); break;
      
    case 2 : dif += logBinom(tryz,dt1->alpha,dt1->beta,dt1->cat,dt1->p,dt1->k) 
	- logBinom(zi,dt1->alpha,dt1->beta,dt1->cat,dt1->p,dt1->k); break;

    case 3 : dif += logPoisson(tryz,dt1->alpha,dt1->beta,dt1->cat,dt1->p,dt1->k) 
	- logPoisson(zi,dt1->alpha,dt1->beta,dt1->cat,dt1->p,dt1->k); break;

    case 4 : dif += logMult(tryz,dt1->alpha,dt1->beta,dt1->cat,dt1->class,dt1->nclass,dt1->p,dt1->c,dt1->k)
	- logMult(zi,dt1->alpha,dt1->beta,dt1->cat,dt1->class,dt1->nclass,dt1->p,dt1->c,dt1->k); break;
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

    case 4 : dif += logMult(tryz,dt2->alpha,dt2->beta,dt2->cat,dt2->class,dt2->nclass,dt2->p,dt2->c,dt2->k)
	- logMult(zi,dt2->alpha,dt2->beta,dt2->cat,dt2->class,dt2->nclass,dt2->p,dt2->c,dt2->k); break;
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
      
    case 4 : dif += logMult(tryz,dt3->alpha,dt3->beta,dt3->cat,dt3->class,dt3->nclass,dt3->p,dt3->c,dt3->k)
	- logMult(zi,dt3->alpha,dt3->beta,dt3->cat,dt3->class,dt3->nclass,dt3->p,dt3->c,dt3->k); break;
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
      
    case 4 : dif += logMult(tryz,dt4->alpha,dt4->beta,dt4->cat,dt4->class,dt4->nclass,dt4->p,dt4->c,dt4->k)
	- logMult(zi,dt4->alpha,dt4->beta,dt4->cat,dt4->class,dt4->nclass,dt4->p,dt4->c,dt4->k); break;
    }
  }

  dif = dif - 0.5*F77_CALL(ddot)(&k,tryz,&incx,tryz,&incy) + 0.5*F77_CALL(ddot)(&k,zi,&incx,zi,&incy);
  
  if((dif > 0) | (runif(0,1) < exp(dif))){
    dvcopy(zi, tryz, k);
    *accept = *accept + 1;
  }

  Free(tryz); 
}

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
	     int *ty0,int *p0,int *c0,double *a0,double *b0,double *con0,int *cat0,int *class0,int *nclass0,
	     double *sigma0,
	     int *ty1,int *p1,int *c1,double *a1,double *b1,double *con1,int *cat1,int *class1,int *nclass1,
	     double *sigma1,
	     int *ty2,int *p2,int *c2,double *a2,double *b2,double *con2,int *cat2,int *class2,int *nclass2,
	     double *sigma2,
	     int *ty3,int *p3,int *c3,double *a3,double *b3,double *con3,int *cat3,int *class3,int *nclass3,
	     double *sigma3, int *accept){

  double **newz, **tempz; /* con22[2][2] = {0,0,0,0}; */
  double **conList[4];
  int **catList[4]; /* cat22[2][2] = {0,0,0,0}; */ /* con22 and cat22 is used for preventing segmentation fault */
  int i, j, ID;
  int conp[4] = {1,1,1,1}; 
  int conn[4] = {1,1,1,1};
  int catp[4] = {1,1,1,1};
  int catn[4] = {1,1,1,1};
  dataType *dt[4];

  /* decide the size of allocated matrix */
  if((*ty0) == 1){
    conp[0] = *p0;
    conn[0] = *n;
  }else{
    catp[0] = *p0;
    catn[0] = *n;
  }

  if((*ty1) == 1){
    conp[1] = *p1;
    conn[1] = *n;
  }else{
    catp[1] = *p1;
    catn[1] = *n;
  }

  if((*ty2) == 1){
    conp[2] = *p2;
    conn[2] = *n;
  }else{
    catp[2] = *p2;
    catn[2] = *n;
  }

  if((*ty3) == 1){
    conp[3] = *p3;
    conn[3] = *n;
  }else{
    catp[3] = *p3;
    catn[3] = *n;
  }

  /*  nrow = (*n)*(*draw); */
  newz = drowm(*n, *k);
  tempz = drowm(*n, *k);

  for(i=0; i<4; i++){
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
  for(j=0; j<(*k); j++){
    for(i=0; i<(*n); i++){
      tempz[i][j] = lastZ[ID]; 
      ID = ID + 1;
    }
  }

  if((*ndt) > 0){
    if((*ty0) == 1){  /* 1 is normal data */
      /*   conList[0] = drowm(*n,*p0); */
      dvtom(conList[0],con0,*n,*p0); /* transform response data to n by p0 matrix */
      /*  catList[0] = irowm(2,2); */ /* originally set to cat22, causing  assignment from incompatible pointer type   */
    }else{             /* 2,3,4 is integer data */
      /*   catList[0] = irowm(*n,*p0); */
      ivtom(catList[0],cat0,*n,*p0);  /* transform response data to n by p0 matrix */
      /*  conList[0] = drowm(2,2); */  /* originally set to con22, causing  sssignment from incompatible pointer type   */
    }   /*conList[0][0] and catList[0][0] is px1 vector; dt[0]->con or dt[0]->cat will be replaced in MCMC */
    fillData(dt[0],a0,b0,conList[0][0],catList[0][0],class0,nclass0,sigma0,p0,k,c0,ty0); 
  }
  /* for the second data set */
  if((*ndt) > 1){
    if((*ty1) == 1){
      /*   conList[1] = drowm(*n,*p1); */
      dvtom(conList[1],con1,*n,*p1);
      /*  catList[1] = irowm(2,2);    */ 
    }else{
      /*  catList[1] = irowm(*n,*p1); */
      ivtom(catList[1],cat1,*n,*p1);
      /*  conList[1] = drowm(2,2); */
    }
    fillData(dt[1],a1,b1,conList[1][0],catList[1][0],class1,nclass1,sigma1,p1,k,c1,ty1);
  }
  /* for the third data set */  
  if((*ndt) > 2){
    if((*ty2) == 1){
      /* conList[2] = drowm(*n,*p2); */
      dvtom(conList[2],con2,*n,*p2);
      /* catList[2] = irowm(2,2);    */
    }else{
      /*  catList[2] = irowm(*n,*p2); */
      ivtom(catList[2],cat2,*n,*p2);
      /* conList[2] = drowm(2,2); */
    }
    fillData(dt[2],a2,b2,conList[2][0],catList[2][0],class2,nclass2,sigma2,p2,k,c2,ty2);
  }
  /* for the four data set */
  if((*ndt) > 3){
    if((*ty3) == 1){
      /*  conList[3] = drowm(*n,*p3); */
      dvtom(conList[3],con3,*n,*p3);
      /*  catList[3] = irowm(2,2); */   
    }else{
      /*  catList[3] = irowm(*n,*p3); */
      ivtom(catList[3],cat3,*n,*p3);
      /* conList[3] = drowm(2,2); */
    }
    fillData(dt[3],a3,b3,conList[3][0],catList[3][0],class3,nclass3,sigma3,p3,k,c3,ty3);
  }

  /* Important to use GetRNGstate */
  GetRNGstate();
  /* MCMC sampling */
  /* If only one data type */
  if(*ndt == 1){
    /* burn.in + the first draw*/
    for(i=0; i<(*burnin); i++){
      for(j=0; j<(*n); j++){
	if(*(dt[0]->type) == 1){
	  dt[0]->con = conList[0][j];
	}else{
	  dt[0]->cat = catList[0][j];
	}
	metroMix(tempz[j],dt[0],dt[1],dt[2],dt[3],sdev,ndt,&accept[j]);
      }
    }

    /*   printf("good after MCMC burning \n"); */
    /* sampling */
    for(i=0; i<(*draw); i++){
      for(j=0; j<(*n); j++){
	if(*(dt[0]->type) == 1){
	  dt[0]->con = conList[0][j];
	}else{
	  dt[0]->cat = catList[0][j];
	}	
	metroMix(tempz[j],dt[0],dt[1],dt[2],dt[3],sdev,ndt,&accept[j]); 
      }
      dmadd(newz,tempz,*n,*k);
    }
  }

  if(*ndt == 2){
    /* burn.in */
    for(i=0; i<(*burnin); i++){
      for(j=0; j<(*n); j++){
	if(*(dt[0]->type) == 1){
	  dt[0]->con = conList[0][j];
	}else{
	  dt[0]->cat = catList[0][j];
	}
	
	if(*(dt[1]->type) == 1){
	  dt[1]->con = conList[1][j];
	}else{
	  dt[1]->cat = catList[1][j];
	}	
	metroMix(tempz[j],dt[0],dt[1],dt[2],dt[3],sdev,ndt,&accept[j]);
      }
    }
    /* sampling */
    for(i=0; i<(*draw); i++){
      for(j=0; j<(*n); j++){
	if(*(dt[0]->type) == 1){
	  dt[0]->con = conList[0][j];
	}else{
	  dt[0]->cat = catList[0][j];
	}

	if(*(dt[1]->type) == 1){
	  dt[1]->con = conList[1][j];
	}else{
	  dt[1]->cat = catList[1][j];
	}	
	metroMix(tempz[j],dt[0],dt[1],dt[2],dt[3],sdev,ndt,&accept[j]);  
      }
      dmadd(newz,tempz,*n,*k);
    }
  }

  if(*ndt == 3){
    /* burn.in */
    for(i=0; i<(*burnin); i++){
      for(j=0; j<(*n); j++){
	if(*(dt[0]->type) == 1){
	  dt[0]->con = conList[0][j];
	}else{
	  dt[0]->cat = catList[0][j];
	}
	
	if(*(dt[1]->type) == 1){
	  dt[1]->con = conList[1][j];
	}else{
	  dt[1]->cat = catList[1][j];
	}

	if(*(dt[2]->type) == 1){
	  dt[2]->con = conList[2][j];
	}else{
	  dt[2]->cat = catList[2][j];
	}
	
	metroMix(tempz[j],dt[0],dt[1],dt[2],dt[3],sdev,ndt,&accept[j]);
      }
    }
    /* sampling */
    for(i=0; i<(*draw); i++){
      for(j=0; j<(*n); j++){
	if(*(dt[0]->type) == 1){
	  dt[0]->con = conList[0][j];
	}else{
	  dt[0]->cat = catList[0][j];
	}

	if(*(dt[1]->type) == 1){
	  dt[1]->con = conList[1][j];
	}else{
	  dt[1]->cat = catList[1][j];
	}

	if(*(dt[2]->type) == 1){
	  dt[2]->con = conList[2][j];
	}else{
	  dt[2]->cat = catList[2][j];
	}
	
	metroMix(tempz[j],dt[0],dt[1],dt[2],dt[3],sdev,ndt,&accept[j]);  
      }
      dmadd(newz,tempz,*n,*k);
    }
  }

  if(*ndt == 4){
    /* burn.in */
    for(i=0; i<(*burnin); i++){
      for(j=0; j<(*n); j++){
	if(*(dt[0]->type) == 1){
	  dt[0]->con = conList[0][j];
	}else{
	  dt[0]->cat = catList[0][j];
	}
	
	if(*(dt[1]->type) == 1){
	  dt[1]->con = conList[1][j];
	}else{
	  dt[1]->cat = catList[1][j];
	}

	if(*(dt[2]->type) == 1){
	  dt[2]->con = conList[2][j];
	}else{
	  dt[2]->cat = catList[2][j];
	}
	
	if(*(dt[3]->type) == 1){
	  dt[3]->con = conList[3][j];
	}else{
	  dt[3]->cat = catList[3][j];
	}

	metroMix(tempz[j],dt[0],dt[1],dt[2],dt[3],sdev,ndt,&accept[j]);
      }
    }
    /* sampling */
    for(i=0; i<(*draw); i++){
      for(j=0; j<(*n); j++){
	if(*(dt[0]->type) == 1){
	  dt[0]->con = conList[0][j];
	}else{
	  dt[0]->cat = catList[0][j];
	}

	if(*(dt[1]->type) == 1){
	  dt[1]->con = conList[1][j];
	}else{
	  dt[1]->cat = catList[1][j];
	}

	if(*(dt[2]->type) == 1){
	  dt[2]->con = conList[2][j];
	}else{
	  dt[2]->cat = catList[2][j];
	}
	
	if(*(dt[3]->type) == 1){
	  dt[3]->con = conList[3][j];
	}else{
	  dt[3]->cat = catList[3][j];
	}

	metroMix(tempz[j],dt[0],dt[1],dt[2],dt[3],sdev,ndt,&accept[j]);  
      }
      dmadd(newz,tempz,*n,*k);
    }
  }

  PutRNGstate();

  dmscale(newz,*n,*k,1.0/(*draw)); /* average of Z */
  dmtov(meanZ,newz,*n,*k);
  dmtov(lastZ,tempz,*n,*k);
  dmfree(newz, *n);
  dmfree(tempz, *n);
  /*
  if((*ndt)>0){
    if((*ty0)==1){
      dmfree(conList[0],(*n));
      imfree(catList[0],2);
    }else{
      imfree(catList[0],(*n));
      dmfree(conList[0],2);
    }
  }

  if((*ndt)>1){
    if((*ty1)==1){
      dmfree(conList[1],(*n));
      imfree(catList[1],2); 
    }else{
      imfree(catList[1],(*n));
      dmfree(conList[1],2);
    }
  }


  if((*ndt)>2){
    if((*ty2)==1){
      dmfree(conList[2],(*n));
      imfree(catList[2],2);
    }else{
      imfree(catList[2],(*n));
      dmfree(conList[2],2);
    }
  }

  if((*ndt)>3){
    if((*ty3)==1){
      dmfree(conList[3],(*n));
      imfree(catList[3],2);
    }else{
      imfree(catList[3],(*n));
      dmfree(conList[3],2);
    }
  }
  */
  for(i=0; i<4; i++){
    dmfree(conList[i],conn[i]);
    imfree(catList[i],catn[i]);
    free(dt[i]);      
  }

 /*  imfree(cat22,2); */
 /*  dmfree(con22,2); */

}


/*
meanZ, lastZ are the output
initZ is the input;
*/
void mcmcBinom(double *meanZ,double *lastZ,double *alpha,double *beta,int *X, double *sigma, double *initz,
		int *n,int *p,int *k,int *burnin,int *draw){
  double **newz, **tempz;
  int **xm;  /*xm is the matrix format of X */
  int i, j, ID;

  /* nrow = (*n)*(*draw); */
  newz = drowm(*n, *k);
  tempz = drowm(*n, *k);
  xm = irowm(*n,*p);
  ivtom(xm,X,(*n),(*p)); /*copy X to xm */

  ID = 0;
  for(j=0; j<(*k); j++){
    for(i=0; i<(*n); i++){
      tempz[i][j] = initz[ID]; 
      /*      tempz[i][j] = 0; */
      ID = ID + 1;
    }
  }

  GetRNGstate();
  /* burn.in + the first draw*/
  for(i=0; i<(*burnin); i++){
    for(j=0; j<(*n); j++){
      metroBinom(tempz[j],alpha,beta,xm[j],p,k,sigma,tempz[j]);
    }
    /*    if((i%1000) == 0){
      printf("- %d -\n",i);
      } */
  }
  /*  printf(" - good after burnin -\n"); */
  /* sampling */
  for(i=0; i<(*draw); i++){
    for(j=0; j<(*n); j++){
      metroBinom(tempz[j],alpha,beta,xm[j],p,k,sigma,tempz[j]);
    }
    dmadd(newz,tempz,*n,*k);
    /*    if((i%100) == 0){
      printf("- %d -\n",i);
      } */
  }
  PutRNGstate();
  /*  printf(" - good after sampling -\n"); */
  dmscale(newz,*n,*k,1.0/(*draw));
  dmtov(meanZ,newz,*n,*k);
  dmtov(lastZ,tempz,*n,*k);
  imfree(xm, *n);
  dmfree(newz, *n);
  dmfree(tempz, *n);
}

/*
meanZ, lastZ are the output
initZ is the input;

 */
void mcmcPoisson(double *meanZ,double *lastZ,double *alpha,double *beta,int *X, double *sigma, double *initz,
		int *n,int *p,int *k,int *burnin,int *draw){
  double **newz, **tempz;
  int **xm;  /*xm is the matrix format of X */
  int i, j, ID;

  /* nrow = (*n)*(*draw); */
  newz = drowm(*n, *k);
  tempz = drowm(*n, *k);
  xm = irowm(*n,*p);
  ivtom(xm,X,(*n),(*p)); /*copy X to xm */

  ID = 0;
  for(j=0; j<(*k); j++){
    for(i=0; i<(*n); i++){
      tempz[i][j] = initz[ID];  
      /* tempz[i][j] = 0; */
      ID = ID + 1;
    }
  }

  GetRNGstate();
  /* burn.in + the first draw*/
  for(i=0; i<(*burnin); i++){
    for(j=0; j<(*n); j++){
      metroPoisson(tempz[j],alpha,beta,xm[j],p,k,sigma,tempz[j]);
    }
    /*    if((i%1000) == 0){
      printf("- %d -\n",i);
      } */
  }
  /*  printf(" - good after burnin -\n"); */
  /* sampling */
  for(i=0; i<(*draw); i++){
    for(j=0; j<(*n); j++){
      metroPoisson(tempz[j],alpha,beta,xm[j],p,k,sigma,tempz[j]);
    }
    dmadd(newz,tempz,*n,*k);
    /*    if((i%100) == 0){
      printf("- %d -\n",i);
      } */
  }
  PutRNGstate();

  /*  printf(" - good after sampling -\n"); */
  dmscale(newz,*n,*k,1.0/(*draw));
  dmtov(meanZ,newz,*n,*k);
  dmtov(lastZ,tempz,*n,*k);
  imfree(xm, *n);
  dmfree(newz, *n);
  dmfree(tempz, *n);
}

/*sigma2 is the estimated variances */

void mcmcNormal(double *meanZ,double *lastZ,double *mu,double *beta,double *sigma2,double *Y, double *sigma, 
		double *initz, int *n,int *p,int *k,int *burnin,int *draw){
  double **newz, **tempz;
  double **ym;  /*ym is the matrix format of Y */
  int i, j; /*, nrow; */

  /* nrow = (*n)*(*draw); */
  newz = drowm(*n, *k);
  tempz = drowm(*n, *k);
  ym = drowm(*n,*p);
  dvtom(ym,Y,(*n),(*p)); /*copy Y to ym */

  for(j=0; j<(*k); j++){
    for(i=0; i<(*n); i++){
      /*      tempz[i][j] = 0; */
      tempz[i][j] = initz[(*n)*j+i];  
    }
  }

  GetRNGstate();
  /* burn.in + the first draw*/
  for(i=0; i<(*burnin); i++){
    for(j=0; j<(*n); j++){
      metroNormal(tempz[j],mu,beta,sigma2,ym[j],p,k,sigma,tempz[j]);
    }
    /*    if((i%1000) == 0){
      printf("- %d -\n",i);
      } */
  }
  /*  printf(" - good after burnin -\n"); */
  /* sampling */
  for(i=0; i<(*draw); i++){
    for(j=0; j<(*n); j++){
      metroNormal(tempz[j],mu,beta,sigma2,ym[j],p,k,sigma,tempz[j]);
    }
    dmadd(newz,tempz,*n,*k);
    /*    if((i%100) == 0){
      printf("- %d -\n",i);
      } */
  }
  PutRNGstate();

  /*  printf(" - good after sampling -\n"); */
  dmscale(newz,*n,*k,1.0/(*draw));
  dmtov(meanZ,newz,*n,*k);
  dmtov(lastZ,tempz,*n,*k);
  dmfree(ym, *n);
  dmfree(newz, *n);
  dmfree(tempz, *n);
}

/* if the data type is multinomial, beta = p  x (k x C) is the format in R; 
   beta must be transposed using t(beta) before passing to this function */
void mcmcMult(double *meanZ,double *lastZ,double *alpha,double *beta,int *X, int *class,int *nclass,
	      double *sigma,double *initz,int *n,int *p, int *c, int *k,int *burnin,int *draw){
  double **newz, **tempz;
  int **xm;
  int i, j, ID;
  /*  double ** mAlpha; */

  /* nrow = (*n)*(*draw); */
  newz = drowm(*n, *k);
  tempz = drowm(*n, *k);
  xm = irowm(*n,*p);
  ivtom(xm,X,(*n),(*p)); /*copy X to xm */
  
  /*  mAlpha = drowm(*p, *c);
    dvtom(mAlpha,alpha,*p,*c);  */

  ID = 0;
  for(j=0; j<(*k); j++){
    for(i=0; i<(*n); i++){
      tempz[i][j] = initz[ID];  
      ID = ID + 1;
    }
  }

  GetRNGstate();
  /*  printf(" - good before burnin -\n");   */
  /* burn.in + the first draw*/
  for(i=0; i<(*burnin); i++){
    for(j=0; j<(*n); j++){
      metroMult(tempz[j],alpha,beta,xm[j],class,nclass,p,c,k,sigma,tempz[j]);
    }
    /*    if((i%1000) == 0){
      printf("- %d -\n",i);
      } */
  }

  /* sampling */
  for(i=0; i<(*draw); i++){
    for(j=0; j<(*n); j++){
      metroMult(tempz[j],alpha,beta,xm[j],class,nclass,p,c,k,sigma,tempz[j]);
    }
    dmadd(newz,tempz,*n,*k);
    /*    if((i%100) == 0){
      printf("- %d -\n",i);
      } */
  }
  PutRNGstate();

  /*  printf(" - good after sampling -\n"); */
  dmscale(newz,*n,*k,1.0/(*draw));
  dmtov(meanZ,newz,*n,*k);
  dmtov(lastZ,tempz,*n,*k);
  dmfree(newz, *n);
  dmfree(tempz, *n);
  imfree(xm, *n);
  /*  dmfree(mAlpha, *p); */
  /*  imfree(mUpc, *p); */
}

