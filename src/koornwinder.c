/* This file is part of SQT, a library for Singular Quadrature on Triangles
 *
 * Copyright (C) 2020 Michael Carley
 *
 * SQT is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  SQT is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SQT.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <glib.h>

#include <blaswrap.h>

#include <sqt.h>

#include "sqt-private.h"

gint SQT_FUNCTION_NAME(sqt_koornwinder_nm)(gint N, SQT_REAL u, SQT_REAL v,
					   gint str, gint imax,
					   SQT_REAL *Knm)

/*
  Definition of Koornwinder polynomials from Greengard et al, Fast
  multipole methods for the evaluation of layer potentials with
  locally-corrected quadratures arXiv:2006.02545v1.

  Jacobi and Legendre recursions from DLMF

  indexing idx_{nm} = n*(n+1)/2 + m
*/
  
{
  gint nn, n, m, idx, idxm1, idxp1, idxP ;
  SQT_REAL x, Jnm1, Jn, tmp ;

  /*initialize (1-v)^m P_m for n = m*/
  x = (2.0*u + v - 1.0)/(1.0 - v) ;
  n = 0 ; m = 0 ;
  idx = n*(n+1)/2 + m ;
  Knm[str*idx] = 1.0 ;
  n = 1 ; m = 1 ;
  idxp1 = n*(n+1)/2 + m ;
  Knm[str*idxp1] = (1.0 - v)*x ;
  
  for ( m = 1 ; m <= N ; m ++ ) {
    n = m ; 
    idxm1 = idx ; idx = idxp1 ; 
    idxp1 = (n+1)*(n+1+1)/2 + m + 1 ;
    Knm[str*idxp1] = (1.0-v)*((2.0*n+1)/(n+1)*x*Knm[str*idx] -
			      ((1.0-v)*n)/(n+1)*Knm[str*idxm1]) ;
  }
  
  x = 1.0 - 2.0*v ;
  for ( m = 0 ; (m < N) ; m ++ ) {
    idxP = m*(m+1)/2 + m ;
    Jnm1 = 1.0 ;
    nn = m ;
    idx = nn*(nn+1)/2 + m ;
    Knm[str*idx] = Jnm1*Knm[str*idxP]*sqrt(2.0*(1+2*m)*(m+1)) ;

    n = 0 ; 
    Jn   = (SQT_REAL)(n+m+1)*(2*n+2*m+3)/(n+1)/(n+2*m+2)*x -
      (SQT_REAL)(2*m+1)*(2*m+1)*(n+m+1)/(n+1)/(n+2*m+2)/(2*n+2*m+1) ;

    nn = m+1 ;
    idx = nn*(nn+1)/2 + m ;
    Knm[str*idx] = Jn*Knm[str*idxP]*sqrt((SQT_REAL)(nn+1)/(m+1)) ;

    for ( n = 1 ; (n <= N-m-1) ; n ++ ) {
      tmp = ((SQT_REAL)(n+m+1)*(2*n+2*m+3)/(n+1)/(n+2*m+2)*x -
	     (SQT_REAL)(2*m+1)*(2*m+1)*(n+m+1)/(n+1)/(n+2*m+2)/(2*n+2*m+1))*Jn -
	(SQT_REAL)n*(n+2*m+1)*(2*n+2*m+3)/(n+1)/(n+2*m+2)/(2*n+2*m+1)*Jnm1 ;
      Jnm1 = Jn ; Jn = tmp ;

      nn = n+m+1 ;
      idx = nn*(nn+1)/2 + m ;
      Knm[str*idx] = Jn*Knm[str*idxP]*sqrt((SQT_REAL)(nn+1)/(m+1)) ;
    }
  }  

  m = N ; nn = N ;
  idx = nn*(nn+1)/2 + m ;
  Knm[str*idx] *= sqrt(2.0*(1+2*m)*(nn+1)) ;

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_koornwinder_interp_matrix)(SQT_REAL *q, gint nq,
						      SQT_REAL *A)

{
  gint i, N ;

  N = 0 ;
  while ( N*(N+1)/2 < nq ) N ++ ;
  
  for ( i = 0 ; i < nq ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_koornwinder_nm)(N, q[i*3+0], q[i*3+1], nq, nq,
					  &(A[i])) ;
#ifdef SQT_SINGLE_PRECISION
    g_assert_not_reached() ;
#else /*SQT_SINGLE_PRECISION*/
    blaswrap_dscal(nq, (q[i*3+2]), &(A[i]), nq) ;
#endif /*SQT_SINGLE_PRECISION*/
  }

  /*return order of Knm required for interpolation after application
    of matrix*/
  return N ;
}
