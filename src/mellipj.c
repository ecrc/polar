
/**
 ZOLOPD-SVD
 *
 * (C) Copyright 2016 King Abdullah University of Science and Technology
 * Authors:
 * Dalal Sukkari (dalal.sukkari@kaust.edu.sa)
 * David Keyes (david.keyes@kaust.edu.sa)
 * Hatem Ltaief (hatem.ltaief@kaust.edu.sa)
 *  
 * Redistribution  and  use  in  source and binary forms, with or without
 * modification,  are  permitted  provided  that the following conditions
 * are met:
 * 
 * Redistributions  of  source  code  must  retain  the above copyright
 * notice,  this  list  of  conditions  and  the  following  disclaimer.
 * Redistributions  in  binary  form must reproduce the above copyright
 * notice,  this list of conditions and the following disclaimer in the
 * documentation  and/or other materials provided with the distribution.
 * Neither  the  name of the King Abdullah University of Science and
 * Technology nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior
 * written permission.
 *
 *
 THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
 LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
function [sn,cn,dn] = mellipj(u,alpha,tol)
mELLIPJ Jacobi elliptic functions. MATLAB's built-in code, modified for improved accuracy. 
   [SN,CN,DN] = ELLIPJ(U,M) returns the values of the Jacobi elliptic 
   functions Sn, Cn and Dn, evaluated for corresponding elements of 
   argument U and parameter M.  U and M must be arrays of the same 
   size or either can be scalar.  As currently implemented, M is 
   limited to 0 <= M <= 1. 

   [SN,CN,DN] = ELLIPJ(U,M,TOL) computes the elliptic functions to
   the accuracy TOL instead of the default TOL = EPS.  
*/
#include "polar.h"

int mellipj(double u, double alpha, double *sn, double *cn, double *dn, double *work)
{

    double eps, tol;
    int chunk = 10, i, n, in;
     
    double *a    = work;
    double *b    = a + chunk + 1;
    double *c    = b + chunk + 1;

    double m1, m, mmax;

    memset(a, 0, (chunk + 1)* sizeof(double));
    memset(b, 0, (chunk + 1)* sizeof(double));
    memset(c, 0, (chunk + 1)* sizeof(double));

    u=creal(u);
    alpha=creal(alpha);

    m  = sin(alpha) * sin(alpha);
    m1 = cos(alpha) * cos(alpha);

    eps = LAPACKE_dlamch_work('e');
    tol = eps;

    mmax = 1;

    *cn=u;
    *sn = cn[0]; *dn = sn[0];


    if ((m < 0) || (m > 1)){ 
	fprintf(stderr, "error(ellipj:MOutOfRange)") ;
	return -1;
    }

    // pre-allocate space and augment if needed
    chunk = 10;
    c[0] = creal(sin(alpha));
    b[0] = creal(cos(alpha));
    a[0] = c[0];

    a[0]=1;

    n = 0.;
    i = 1;
    in = 0;
    while ((fabs(c[i-1]) > tol) && i<1000 ){
        i = i + 1;
        if (i > chunk+1){

           a = work;
           b = a + (2*chunk + 1);
           c = b + (2*chunk + 1);
           memset(a, 0, (2*chunk + 1)* sizeof(double));
           memset(b, 0, (2*chunk + 1)* sizeof(double));
           memset(c, 0, (2*chunk + 1)* sizeof(double));
        }
        a[i-1] = 0.5 * (a[i-2] + b[i-2]);
        b[i-1] = sqrt(a[i-2] * b[i-2]);
        c[i-1] = 0.5 * (a[i-2] - b[i-2]);
        if ((fabs(c[i-1]) <= tol) && (fabs(c[i-2]) > tol)) {in = 1;}
        if (in != 0){
           n = i-1;
        }
    }
    double *phin = c + (2*chunk + 1);
    memset(phin, 0, i * sizeof(double));
    phin[0]=u;

    phin[i-1] = cpow(2,n)*a[i-1]*u;
    while (i > 1){
        i = i - 1;
        if (n >= i) {in = 1;}
        phin[i-1] = phin[i];
        if (in != 0){
           phin[i-1] = 0.5 * (asin(c[i]*sin(creal(phin[i]))/a[i]) + phin[i]);
        }
    }
    *sn = sin(creal(phin[0]));
    *cn = cos(creal(phin[0]));
    *dn = sqrt(1 - m * sn[0]*sn[0]);

return 0;
}

int choosem(double con, int *m){  
    //choose the Zolotarev degree according to the condnum
    if (con < 1.001)         {*m=2;  }        
    else if (con <= 1.01)    {*m = 3;}            
    else if (con <= 1.1)     {*m = 4;}            
    else if (con <= 1.2)     {*m = 5;}            
    else if (con <= 1.5)     {*m = 6;}
    else if (con <= 2)       {*m = 8;}  // one-step convergence till here
    else if (con < 6.5)      {*m = 2;}
    else if (con < 180)      {*m = 3;}       
    else if (con < 1.5*1e4) {*m = 4;}
    else if (con < 2*1e6)   {*m = 5;}
    else if (con < 1*1e9)   {*m = 6;}
    else if (con < 3*1e12)  {*m = 7;}
    else                     {*m = 8;}
    
return 0;
}

