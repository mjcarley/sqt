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

static gint koornwinder_deriv_recursion(gint N, gint m,
					SQT_REAL Pm, SQT_REAL dPm,
					SQT_REAL vpm,
					SQT_REAL u, SQT_REAL v,
					SQT_REAL w,
					SQT_REAL *K, gint kstr,
					SQT_REAL *Ku, gint ustr,
					SQT_REAL *Kv, gint vstr,
					gint nst)

{
  SQT_REAL cnm, tmp, Jnm1, Jnm, dJnm, A, B, C ;
  gint n, idx ;

  Jnm1 = 1.0 ; dJnm = 0.0 ;
  n = m ;
  idx = n*(n+1)/2 + m ;
  cnm = sqrt(2.0*(1+2*m)*(m+1)) ;

  if ( idx < nst ) {
    K [kstr*idx] = Jnm1* Pm*vpm*cnm ;
    Ku[ustr*idx] = Jnm1*dPm*vpm*cnm*2.0/(1.0-v) ;
    Kv[vstr*idx] =
      Ku[ustr*idx]*u/(1.0-v) - cnm*vpm*Pm*(Jnm1*m/(1.0-v) + dJnm*2.0) ;
  }
  if ( m+1 > N ) return 0 ;
  
  A = 0.5*(2*m+3) ; B = -0.5*(2*m+1) ;

  Jnm = A*w + B ; dJnm = A ;
  n = m + 1 ;
  cnm = sqrt(2.0*(1+2*m)*(n+1)) ;
  idx = n*(n+1)/2 + m ;
  if ( idx < nst ) {
    K [kstr*idx] = Jnm* Pm*vpm*cnm ;
    Ku[ustr*idx] = Jnm*dPm*vpm*cnm*2.0/(1.0-v) ;
    Kv[vstr*idx] =
      Ku[ustr*idx]*u/(1.0-v) - cnm*vpm*Pm*(Jnm*m/(1.0-v) + dJnm*2.0) ;
  }

  if ( N*(N+1)/2 + m >= nst ) N -- ;

  for ( n = m+2 ; n <= N ; n ++ ) {
    A =  (SQT_REAL)(n-1+1)*(2*n-2+3)/(n-1-m+1)/(n-1+m+2) ;
    B = -(SQT_REAL)(2*m+1)*(2*m+1)*(n-1+1)/(n-1-m+1)/(n-1+m+2)/(2*n-2+1) ;
    C =  (SQT_REAL)(n-1-m)*(n-1+m+1)*(2*n-2+3)/(n-1-m+1)/(n-1+m+2)/(2*n-2+1) ;
    
    tmp = Jnm ;
    Jnm = (A*w + B)*Jnm - C*Jnm1 ;
    Jnm1 = tmp ;

    dJnm  = 2.0*(n-m)*(n+m+1)*Jnm1 - (n-m)*((2*n+1)*w + 2*m+1)*Jnm ;
    dJnm /= (1.0-w*w)*(2*n+1) ;

    cnm = sqrt(2.0*(1+2*m)*(n+1)) ;
    idx = n*(n+1)/2 + m ;
    K [kstr*idx] = Jnm* Pm*vpm*cnm ;
    Ku[ustr*idx] = Jnm*dPm*vpm*cnm*2.0/(1.0-v) ;
    Kv[vstr*idx] =
      Ku[ustr*idx]*u/(1.0-v) - cnm*vpm*Pm*(Jnm*m/(1.0-v) + dJnm*2.0) ;
      /* -cnm*vpm*(Jnm*(Pm*m - 2.0*dPm*u/(1.0-v))/(1.0-v) + Pm*dJnm*2.0) ; */
  }

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_koornwinder_deriv_nm)(gint N, SQT_REAL u, SQT_REAL v,
						 SQT_REAL *K , gint kstr,
						 SQT_REAL *Ku, gint ustr,
						 SQT_REAL *Kv, gint vstr,
						 gint nst)

/*
  Definition of Koornwinder polynomials from Greengard et al, Fast
  multipole methods for the evaluation of layer potentials with
  locally-corrected quadratures arXiv:2006.02545v1.

  Jacobi and Legendre recursions from DL2MF

  indexing idx_{nm} = n*(n+1)/2 + m
*/
  
{
  gint m ;
  SQT_REAL x, tmp, w, Pm, Pmm1, dPm, vpm ;
  
  w = 1.0 - 2.0*v ;
  x = (2.0*u + v - 1.0)/(1.0 - v) ;

  m = 0 ;
  Pm = 1.0 ; dPm = 0.0 ; vpm = 1.0 ; 
  koornwinder_deriv_recursion(N, m, Pm, dPm, vpm, u, v, w,
			      K, kstr, Ku, ustr, Kv, vstr, nst) ;

  m = 1 ;
  Pmm1 = Pm ; Pm = x ; dPm = 1.0 ; vpm = 1.0 - v ;

  koornwinder_deriv_recursion(N, m, Pm, dPm, vpm, u, v, w,
			      K, kstr, Ku, ustr, Kv, vstr, nst) ;

  /*recursive generation of Legendre polynomials*/
  for ( m = 2 ; m <= N ; m ++ ) {
    tmp = Pm ;
    Pm = (Pm*x*(2*m-1) - Pmm1*(m-1))/m ;
    Pmm1 = tmp ;
    dPm = (Pmm1 - x*Pm)/(1.0-x*x)*m ;
    vpm *= 1.0 - v ;
    koornwinder_deriv_recursion(N, m, Pm, dPm, vpm, u, v, w,
				K, kstr, Ku, ustr, Kv, vstr, nst) ;
  }
  
  return 0 ;
}

static gint koornwinder_recursion(gint N, gint m, SQT_REAL Pm, SQT_REAL w,
				  SQT_REAL *Knm, gint str, gint nst)

{
  SQT_REAL cnm, tmp, Jnm1, Jnm, A, B, C ;
  gint n, idx ;

  Jnm1 = 1.0 ;
  n = m ;
  idx = n*(n+1)/2 + m ;
  cnm = sqrt(2.0*(1+2*m)*(m+1)) ;

  if ( idx < nst ) Knm[str*idx] = Jnm1*Pm*cnm ;
  if ( m+1 > N ) return 0 ;
  
  A = 0.5*(2*m+3) ; B = -0.5*(2*m+1) ;

  Jnm = A*w + B ;
  n = m + 1 ;
  cnm = sqrt(2.0*(1+2*m)*(n+1)) ;
  idx = n*(n+1)/2 + m ;
  if ( idx < nst ) Knm[str*idx] = Jnm*Pm*cnm ;

  if ( N*(N+1)/2 + m >= nst ) N -- ;

  for ( n = m+2 ; n <= N ; n ++ ) {
    A =  (SQT_REAL)(n-1+1)*(2*n-2+3)/(n-1-m+1)/(n-1+m+2) ;
    B = -(SQT_REAL)(2*m+1)*(2*m+1)*(n-1+1)/(n-1-m+1)/(n-1+m+2)/(2*n-2+1) ;
    C =  (SQT_REAL)(n-1-m)*(n-1+m+1)*(2*n-2+3)/(n-1-m+1)/(n-1+m+2)/(2*n-2+1) ;
    
    tmp = Jnm ;
    Jnm = (A*w + B)*Jnm - C*Jnm1 ;
    Jnm1 = tmp ;
    cnm = sqrt(2.0*(1+2*m)*(n+1)) ;
    idx = n*(n+1)/2 + m ;
    Knm[str*idx] = Jnm*Pm*cnm ;
  }

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_koornwinder_nm)(gint N, SQT_REAL u, SQT_REAL v,
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
  
  w = 1.0 - 2.0*v ;
  x = (2.0*u + v - 1.0)/(1.0 - v) ;

  m = 0 ;
  Pm = 1.0 ;
  koornwinder_recursion(N, m, Pm, w, Knm, str, nst) ;

  m = 1 ;
  Pmm1 = Pm ; Pm = (1.0 - v)*x ;
  koornwinder_recursion(N, m, Pm, w, Knm, str, nst) ;

  /*recursive generation of Legendre polynomials and Jacobi recursion*/
  for ( m = 2 ; m <= N ; m ++ ) {
    tmp = Pm ;
    Pm = (1.0-v)*((2.0*m-1)*x*Pm - (1.0-v)*(m-1)*Pmm1)/m ;
    Pmm1 = tmp ;
    koornwinder_recursion(N, m, Pm, w, Knm, str, nst) ;
  }  
  
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
  N -- ;
  
  for ( i = 0 ; i < nst ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_koornwinder_nm)(N, s[i*sstr], t[i*tstr],
    					  &(A[i]), nst, nst) ;
#ifdef SQT_SINGLE_PRECISION
    g_assert_not_reached() ;
#else /*SQT_SINGLE_PRECISION*/
    blaswrap_dscal(nst, (w[i*wstr]), &(A[i]), nst) ;
#endif /*SQT_SINGLE_PRECISION*/
  }
  
  /*return order of Knm required for interpolation after application
    of matrix

    matrix contains polynomials K_nm, 0 <= n <= N, 0 <= m <= n, 
    n*(n+1)/2 + m < nst
  */

  return N ;
}
