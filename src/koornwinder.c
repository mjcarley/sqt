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

/* #define INDEX_TRACKING */

gint SQT_FUNCTION_NAME(sqt_koornwinder_nm)(gint N, SQT_REAL u, SQT_REAL v,
					   SQT_REAL *Knm, gint str)

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
#ifdef INDEX_TRACKING
  gint imax = 0 ;
#endif /*INDEX_TRACKING*/
  
  /*initialize (1-v)^m P_m for n = m*/
  x = (2.0*u + v - 1.0)/(1.0 - v) ;
  n = 0 ; m = 0 ;
  idx = n*(n+1)/2 + m ;
#ifdef INDEX_TRACKING
  imax = idx ; 
#endif /*INDEX_TRACKING*/
  Knm[str*idx] = 1.0 ;
  n = 1 ; m = 1 ;
  idxp1 = n*(n+1)/2 + m ;
  Knm[str*idxp1] = (1.0 - v)*x ;
#ifdef INDEX_TRACKING
  imax = MAX(imax, idxp1) ; 
#endif /*INDEX_TRACKING*/

  /*recursive generation of Legendre polynomials*/
  for ( m = 1 ; m < N ; m ++ ) {
    n = m ; 
    idxm1 = idx ; idx = idxp1 ; 
    idxp1 = (n+1)*(n+1+1)/2 + m + 1 ;
    Knm[str*idxp1] = (1.0-v)*((2.0*n+1)/(n+1)*x*Knm[str*idx] -
			      ((1.0-v)*n)/(n+1)*Knm[str*idxm1]) ;
#ifdef INDEX_TRACKING
  imax = MAX(imax, idxp1) ; 
#endif /*INDEX_TRACKING*/    
  }
  
  x = 1.0 - 2.0*v ;
  for ( m = 0 ; (m < N) ; m ++ ) {
    /*index of Legendre polynomial P_m*/
    idxP = m*(m+1)/2 + m ;
    Jnm1 = 1.0 ;
    nn = m ;
    idx = nn*(nn+1)/2 + m ;
    Knm[str*idx] = Jnm1*Knm[str*idxP]*sqrt(2.0*(1+2*m)*(m+1)) ;
#ifdef INDEX_TRACKING
  imax = MAX(imax, idx) ; 
#endif /*INDEX_TRACKING*/

    n = 0 ; 
    Jn   = (SQT_REAL)(n+m+1)*(2*n+2*m+3)/(n+1)/(n+2*m+2)*x -
      (SQT_REAL)(2*m+1)*(2*m+1)*(n+m+1)/(n+1)/(n+2*m+2)/(2*n+2*m+1) ;

    nn = m+1 ;
    idx = nn*(nn+1)/2 + m ;
    Knm[str*idx] = Jn*Knm[str*idxP]*sqrt((SQT_REAL)(nn+1)/(m+1)) ;
#ifdef INDEX_TRACKING
  imax = MAX(imax, idx) ; 
#endif /*INDEX_TRACKING*/

    for ( n = 1 ; (n <= N-m-1) ; n ++ ) {
      /*Jacobi polynomial recursion*/
      tmp = ((SQT_REAL)(n+m+1)*(2*n+2*m+3)/(n+1)/(n+2*m+2)*x -
	     (SQT_REAL)(2*m+1)*(2*m+1)*(n+m+1)/(n+1)/(n+2*m+2)/(2*n+2*m+1))*Jn -
	(SQT_REAL)n*(n+2*m+1)*(2*n+2*m+3)/(n+1)/(n+2*m+2)/(2*n+2*m+1)*Jnm1 ;
      Jnm1 = Jn ; Jn = tmp ;

      nn = n+m+1 ;
      idx = nn*(nn+1)/2 + m ;
      Knm[str*idx] = Jn*Knm[str*idxP]*sqrt((SQT_REAL)(nn+1)/(m+1)) ;
#ifdef INDEX_TRACKING
  imax = MAX(imax, idx) ; 
#endif /*INDEX_TRACKING*/
    }
  }  

  m = N ; nn = N ;
  idx = nn*(nn+1)/2 + m ;
  Knm[str*idx] *= sqrt(2.0*(1+2*m)*(nn+1)) ;
#ifdef INDEX_TRACKING
  imax = MAX(imax, idx) ;

  fprintf(stderr, "maximum index: %d\n", imax) ;
#endif /*INDEX_TRACKING*/

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_koornwinder_deriv_nm)(gint N, SQT_REAL u, SQT_REAL v,
						 SQT_REAL *K , gint kstr,
						 SQT_REAL *Ku, gint ustr,
						 SQT_REAL *Kv, gint vstr)

/*
  Evaluation of Koornwinder polynomials and derivatives

  Definition of Koornwinder polynomials from Greengard et al, Fast
  multipole methods for the evaluation of layer potentials with
  locally-corrected quadratures arXiv:2006.02545v1.

  Jacobi and Legendre recursions from DLMF

  indexing idx_{nm} = n*(n+1)/2 + m
*/
  
{
  gint n, m, idx, idxm1, idxp1, idxP ;
  SQT_REAL x, Jnm1, Jn, dJn, tmp, Pm, dPm, cnm, vpm, A, B, C ;

  idx = N*(N+1)/2 ;
  
  /*initialize P_m for n = m*/
  x = (2.0*u + v - 1.0)/(1.0 - v) ;
  n = 0 ; m = 0 ;
  idx = n*(n+1)/2 + m ;
  K[kstr*idx] = 1.0 ; Ku[ustr*idx  ] = 0.0 ; Kv[vstr*idx  ] = 0.0 ; 
  n = 1 ; m = 1 ;
  idxp1 = n*(n+1)/2 + m ;
  K[kstr*idxp1] = x ; Ku[ustr*idxp1] = 1.0 ; 

  /*initialize Legendre polynomials and derivatives*/
  for ( m = 1 ; m < N ; m ++ ) {
    n = m ; 
    idxm1 = idx ; idx = idxp1 ; 
    idxp1 = (n+1)*(n+1+1)/2 + m + 1 ;
    K[kstr*idxp1] = (2.0*m+1)/(m+1)*x*K[kstr*idx] -
      (SQT_REAL)(m)/(m+1)*K[kstr*idxm1] ;

    /*derivative of P_m*/
    Ku[ustr*idxp1] = (K[kstr*idx] - x*K[kstr*idxp1])*(m+1)/(1.0 - x*x) ;
  }
  
  x = 1.0 - 2.0*v ; vpm = 1.0 ;
  for ( m = 0 ; (m <= N) ; m ++ ) {
    idxP = m*(m+1)/2 + m ;
    Pm  = K[kstr*idxP] ; dPm = Ku[ustr*idxP] ;
    Jnm1 = 1.0 ; dJn = 0.0 ;
    n = m ;
    idx = n*(n+1)/2 + m ;
    cnm = sqrt(2.0*(1+2*m)*(n+1)) ;
    K [kstr*idx] = cnm*vpm*Jnm1*Pm ;
    Ku[ustr*idx] = cnm*vpm*Jnm1*dPm*2.0/(1.0-v) ;
    Kv[vstr*idx] =
      -cnm*vpm*(Jnm1*(Pm*m - 2.0*dPm*u/(1.0-v))/(1.0-v) + Pm*dJn*2.0) ;

    A = 0.5*(2*m+3) ; B = -0.5*(2*m+1) ;

    Jn = A*x + B ; dJn = A ;
    
    n = m+1 ;
    cnm = sqrt((2.0)*(2*m+1)*(n+1)) ;
    idx = n*(n+1)/2 + m ;
    K[kstr*idx] = cnm*vpm*Jn*Pm ;
    Ku[ustr*idx] = cnm*vpm*Jn*dPm*2.0/(1.0-v) ;
    Kv[vstr*idx] =
	-cnm*vpm*(Jn*(Pm*m - 2.0*dPm*u/(1.0-v))/(1.0-v) + Pm*dJn*2.0) ;

    for ( n = m+2 ; n <= N ; n ++ ) {
      A =  (SQT_REAL)((n)*(2*n+1))/(n-m)/(n+m+1) ;
      B = -(SQT_REAL)((2*m+1)*(2*m+1)*(n))/(n-m)/(n+m+1)/(2*n-1) ;
      C =  (SQT_REAL)((n-m-1)*(n+m)*(2*n+1))/(n-m)/(n+m+1)/(2*n-1) ;
      
      tmp = (A*x + B)*Jn - C*Jnm1 ;
      
      Jnm1 = Jn ; Jn = tmp ;
      dJn = 2.0*(n-m)*(n+m+1)*Jnm1 - (n-m)*((2*n+1)*x + 2*m+1)*Jn ;
      dJn /= (2*n+1)*(1.0-x*x) ;
      
      idx = n*(n+1)/2 + m ;
      cnm = sqrt(2.0*(1+2*m)*(n+1)) ;
      K[kstr*idx] = cnm*vpm*Jn*Pm ;
      Ku[ustr*idx] = cnm*vpm*Jn*dPm*2.0/(1.0-v) ;
      Kv[vstr*idx] =
	-cnm*vpm*(Jn*(Pm*m - 2.0*dPm*u/(1.0-v))/(1.0-v) + Pm*dJn*2.0) ;
    }

    vpm *= 1.0 - v ;
  }  

  return 0 ;
}

static gint koornwinder_recursion(gint N, gint m, SQT_REAL Pm, SQT_REAL w,
				  SQT_REAL *Knm, gint str, gint nst)

{
  SQT_REAL cnm, tmp, Jnm1, Jn ;
  gint n, idx, nmax ;

  nmax = N-m-1 ;
  nmax = (gint)floor(-m-1.5+0.5*sqrt(1 + 8*(nst -1 - m))) ;

  if ( nmax <= 0 ) return 0 ;
  Jnm1 = 1.0 ;
  n = m ;
  idx = n*(n+1)/2 + m ;
  /* g_assert(idx < nst) ; */
  cnm = sqrt(2.0*(1+2*m)*(m+1)) ;
  Knm[str*idx] = Jnm1*Pm*cnm ;
  
  n = 0 ; 
  Jn = (SQT_REAL)(n+m+1)*(2*n+2*m+3)/(n+1)/(n+2*m+2)*w -
    (SQT_REAL)(2*m+1)*(2*m+1)*(n+m+1)/(n+1)/(n+2*m+2)/(2*n+2*m+1) ;

  if ( nmax < 1 ) return 0 ;
  n = m+1 ;
  idx = n*(n+1)/2 + m ;
  /* g_assert(idx < nst) ; */
  cnm = sqrt(2.0*(1+2*m)*(n+1)) ;
  Knm[str*idx] = Jn*Pm*cnm ;

  for ( (n = 1), (idx = (n+m+1)*(n+m+2)/2 + m) ;
	(n <= nmax) ;
	(n ++), (idx = (n+m+1)*(n+m+2)/2 + m) ) {
    /* g_assert(idx < nst) ; */
    /*Jacobi polynomial recursion*/
    tmp = ((SQT_REAL)(n+m+1)*(2*n+2*m+3)/(n+1)/(n+2*m+2)*w -
	   (SQT_REAL)(2*m+1)*(2*m+1)*(n+m+1)/(n+1)/(n+2*m+2)/(2*n+2*m+1))*Jn -
      (SQT_REAL)n*(n+2*m+1)*(2*n+2*m+3)/(n+1)/(n+2*m+2)/(2*n+2*m+1)*Jnm1 ;
    Jnm1 = Jn ; Jn = tmp ;
    cnm = sqrt(2.0*(1+2*m)*(n+m+2)) ;
    Knm[str*idx] = Jn*Pm*cnm ;
  }

  return 0 ;
}

static gint sqt_koornwinder_tmp(gint N, SQT_REAL u, SQT_REAL v,
				SQT_REAL *Knm, gint str, gint nst)

/*
  Definition of Koornwinder polynomials from Greengard et al, Fast
  multipole methods for the evaluation of layer potentials with
  locally-corrected quadratures arXiv:2006.02545v1.

  Jacobi and Legendre recursions from DLMF

  indexing idx_{nm} = n*(n+1)/2 + m
*/
  
{
  gint m ;
  SQT_REAL x, tmp, w, Pm, Pmm1 ;
  
  /*initialize (1-v)^m P_m for n = m*/
  w = 1.0 - 2.0*v ;
  x = (2.0*u + v - 1.0)/(1.0 - v) ;

  m = 0 ;
  Pm = 1.0 ;
  koornwinder_recursion(N, m, Pm, w, Knm, str, nst) ;

  m = 1 ;
  Pmm1 = Pm ; Pm = (1.0 - v)*x ;
  koornwinder_recursion(N, m, Pm, w, Knm, str, nst) ;

  /*recursive generation of Legendre polynomials*/
  for ( m = 2 ; m < N ; m ++ ) {
    tmp = Pm ;
    Pm = (1.0-v)*((2.0*m-1)/m*x*Pm - ((1.0-v)*(m-1))/m*Pmm1) ;
    Pmm1 = tmp ;
    koornwinder_recursion(N, m, Pm, w, Knm, str, nst) ;
  }  

  /* exit(0) ; */
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_koornwinder_interp_matrix)(SQT_REAL *s, gint sstr,
						      SQT_REAL *t, gint tstr,
						      SQT_REAL *w, gint wstr,
						      gint nst,
						      SQT_REAL *A)

{
  gint i, N ;

  N = 0 ;
  while ( N*(N+1)/2 < nst ) N ++ ;
  
  for ( i = 0 ; i < nst ; i ++ ) {
    /* SQT_FUNCTION_NAME(sqt_koornwinder_nm)(N, s[i*sstr], t[i*tstr], */
    /* 					  &(A[i]), nst) ; */
    sqt_koornwinder_tmp(N, s[i*sstr], t[i*tstr],
			&(A[i]), nst, nst) ;
#ifdef SQT_SINGLE_PRECISION
    g_assert_not_reached() ;
#else /*SQT_SINGLE_PRECISION*/
    blaswrap_dscal(nst, (w[i*wstr]), &(A[i]), nst) ;
#endif /*SQT_SINGLE_PRECISION*/
  }
  
  /*return order of Knm required for interpolation after application
    of matrix*/
  return N ;
}
