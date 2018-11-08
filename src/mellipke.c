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
 **/
#include "polar.h"

int mellipke( double alpha, double *k, double *e)
{
    /* ELLIPKE Complete elliptic integral. Modified from Matlab's built-in code 
     * for improved accuracy
     */
    double m, m1, mm, a0, b0, s0;
    double a1, b1, c1, w1;
    double eps, tol;
    double mmold;
    int i1;

    m = sin(alpha) * sin(alpha);
    m1 = cos(alpha) * cos(alpha);

    eps = LAPACKE_dlamch_work('e');
    tol = eps;
    

    if (m == 0.0) {k[0] = 0.0; *e = k[0]; return 0;}
    if (m < 0 || m > 1){ 
	fprintf(stderr, "error(ellipke:MOutOfRange)") ;
	return -1;
    }

    a0 = 1.;
    b0 = cos(alpha);
    s0 = m;
    i1 = 0.; mm = 1.;

    while ( mm > tol){
        a1 = (a0+b0)/2;
        b1 = sqrt(a0*b0);
        c1 = (a0-b0)/2;
        i1 = i1 + 1;
        w1 = cpow(2,i1)*c1*c1;
        mm = w1;
        s0 = s0 + w1;
        a0 = a1;
        b0 = b1;
    }

    *k = M_PI/(2*a1);
    *e = k[0]*(1-s0/2);
return 0;
}
