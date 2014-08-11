#include <R.h>
#include <Rinternals.h>

SEXP thex(SEXP x, SEXP y, SEXP wts, SEXP hbin, SEXP xhc, SEXP yhc){
  
  
  int n, prcnt=0;
  SEXP htcnt,  htwcnt;

  PROTECT(htcnt = allocVector(INTSXP, lmax*6));
  PROTECT(ytc = allocVector(REALSXP, lmax*6));
  PROTECT(xtc = allocVector(REALSXP, lmax*6));
  prcnt=prcnt+3;
  if(doWeights>0)PROTECT(htwcnt = allocVector(REALSXP, lmax*6));
  else PROTECT(htwcnt = allocVector(NILSXP, 1));
  prcnt++;	 

  memset(INTEGER(htcnt),0,6*lmax*sizeof(int));
  memset(REAL(xhc),0,lmax*sizeof(double));
  hcx=REAL(xhc);
  memset(REAL(xhc),0,lmax*sizeof(double));
  hcy=REAL(yhc);
  if(doWeights>0) memset(REAL(htwcnt),0,6*lmax*sizeof(double));
  
  hci=floor(ll/jinc);
  hcj=fmod(ll,jinc);
  hcy=c4 * hci + ymin;
  if(fmod(hci,2.0)!=0) hcx=c3*(hcj + 0.5) + xmin;
  else hcx=c3*hcj + xmin;
  if(REAL(x)[i]>hcx) tcell = 3;
	    
}
