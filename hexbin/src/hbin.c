#include <R.h>
#include <Rinternals.h>


/* int and double zeros are all bits off */
#define ZEROINT(X,N,I) do{memset(INTEGER(X),0,N*sizeof(int));}while(0)
#define ZERODBL(X,N,I) do{memset(REAL(X),0,N*sizeof(double));}while(0)


/* 
## an R test
x<-rnom(100);y<-rnom(100)
xbins<-10;
jmax <- floor(xbins + 1.5001)
c1 <- 2 * floor((xbins)/sqrt(3) + 1.5001)
imax <- trunc((jmax*c1 -1)/jmax + 1)
dyn.load(hbin.so)
.Call("hbin",x,y,NULL,1.0,xbins,range(x),range(y),as.integer(c(imax,jmax)),100L,FALSE,FALSE)
*/

SEXP hbin(SEXP x, SEXP y, SEXP swts, SEXP shape, 
          SEXP size, SEXP rx, SEXP ry, SEXP bnd, SEXP n,
	  SEXP doCellid){
/*	Copyright 1991
	Version Date:	September 16, 1994
	Programmer:	Dan Carr, Conversion to C, triangulation, 
                        and Rapi Nicholas Lewin-Koh (2010)
	Indexing:	Left to right, bottom to top
			bnd[0] rows, bnd[2] columns
        Input Vars:
			x,y       the values of x and y
			xcm,ycm   vectors for the center of mass of the returned hexagons
			shape     the shape parameter for the hexagons
			cell
			cnt 
	Output:	     
                     cell ids for non empty cells, revised bnd(1)
                     optionally also return cellid(1:n), and wcnt
      Copyright (2004) Nicholas Lewin-Koh and Martin Maechler */


  int nc, nn;
  int i, i1, i2, iinc;
  int j1, j2, jinc;
  int L, ll, lmax, lat, tcell;
  double c1, c2, con1, con2, dist1, fsize;
  double sx, sy, xmin, ymin, xr, yr;
  uint keepID=0, doWeights=0;
  int prcnt=0;
  SEXP ans;
  SEXP cnt, cell, wcnt, cellid,  xcm,  ycm;
	
        

	if(LOGICAL(doCellid)[0]>0) keepID = 1;
	if(length(swts) > 0 || swts != R_NilValue) doWeights = 1;
       /*_______Alloc and protect the necessary result vectors, then set to 0_____________*/
        lmax=INTEGER(bnd)[0]*INTEGER(bnd)[1];
    /* prcnt is the number of protected R vectors */
	PROTECT(cnt = allocVector(INTSXP, lmax));
        prcnt++;
	PROTECT(cell = allocVector(INTSXP, lmax));
        prcnt++;
	PROTECT(xcm = allocVector(REALSXP, lmax));
        prcnt++;
	PROTECT(ycm = allocVector(REALSXP, lmax));
        prcnt++;	 
	if(keepID > 0) PROTECT(cellid = allocVector(INTSXP, INTEGER(n)[0]));
	else PROTECT(cellid = allocVector(NILSXP, 1));
	prcnt++;
	  
	if(doWeights > 0)PROTECT(wcnt = allocVector(REALSXP, lmax));
	else PROTECT(wcnt = allocVector(NILSXP, 1));
	prcnt++;	


	memset(INTEGER(cell),0,lmax*sizeof(int));
	memset(INTEGER(cnt),0,lmax*sizeof(int));
	memset(REAL(xcm),0,lmax*sizeof(double));
	memset(REAL(ycm),0,lmax*sizeof(double));

	if(doWeights>0) memset(REAL(wcnt),0,lmax*sizeof(double));

       /*_______Constants for scaling the data_____________________________*/
	fsize=INTEGER(size)[0];
	nn=INTEGER(n)[0];
	xmin = REAL(rx)[0];
	ymin = REAL(ry)[0];
	xr = REAL(rx)[1]-xmin;
	yr = REAL(ry)[1]-ymin;
        
	c1 = fsize/xr;
	c2 = (fsize*REAL(shape)[0])/(yr*sqrt(3.0));

	jinc= INTEGER(bnd)[1];
	lat=jinc+1;
	iinc= 2*jinc;
	con1 = 0.25;
	con2 = 1.0/3.0;
        
	/*_______Binning loop______________________________________________*/
	
	for(i=0; i<nn; i++){
	  sx = c1 * (REAL(x)[i] - xmin);
	  sy = c2 * (REAL(y)[i] - ymin);
	  j1 = sx+.5;
	  i1 = sy+.5;
	  dist1 = (sx-j1)*(sx-j1)+ 3.0*(sy-i1)*(sy-i1);
	  /* need floor in C for this, same effect as trunc*/
	  if(dist1 < con1) L = i1*iinc + j1 + 1;
	  else if(dist1 > con2) L = floor(sy)*iinc + floor(sx) + lat;
	  else{
	    j2 = sx;
	    i2 = sy;
	    if(dist1 <= ((sx - j2 - 0.5)*(sx - j2 - 0.5)) + 3.0*((sy - i2 - 0.5)*(sy - i2 - 0.5))) L = i1*iinc + j1 + 1;
	    else L=i2*iinc+ j2+lat;
	  }
	  ll=L-1;
	  INTEGER(cnt)[ll]++;
	  if(doWeights > 0) REAL(wcnt)[ll] = REAL(wcnt)[ll] + REAL(swts)[i];
	  if (keepID > 0) INTEGER(cellid)[i] = L;
	  REAL(xcm)[ll] = REAL(xcm)[ll] + (REAL(x)[i]-REAL(xcm)[ll])/INTEGER(cnt)[ll];
	  REAL(ycm)[ll] = REAL(ycm)[ll]+ (REAL(y)[i]-REAL(ycm)[ll])/INTEGER(cnt)[ll];
	}

/*_______Compression of output________________________________________*/

        nc=-1;
        for(L=0;L<lmax;L++){
	  if(INTEGER(cnt)[L] > 0){
	    nc=nc+1;
	    INTEGER(cell)[nc]=L+1;
	    INTEGER(cnt)[nc]=INTEGER(cnt)[L];	    
	    REAL(xcm)[nc]=REAL(xcm)[L];
	    REAL(ycm)[nc]=REAL(ycm)[L];
	  }
	}

       INTEGER(n)[0]=nc+1;
       INTEGER(bnd)[0]=(INTEGER(cell)[nc]-1)/INTEGER(bnd)[1]+1;
       /* Output constructor */
       PROTECT(ans = allocVector(VECSXP, 6));
       prcnt++;
       SET_VECTOR_ELT(ans, 0, n);
       SET_VECTOR_ELT(ans, 1, cell);
       SET_VECTOR_ELT(ans, 2, cnt);
       SET_VECTOR_ELT(ans, 3, wcnt);
       SET_VECTOR_ELT(ans, 4, xcm);
       SET_VECTOR_ELT(ans, 5, ycm);       

       UNPROTECT(prcnt);
       return(ans);
}
