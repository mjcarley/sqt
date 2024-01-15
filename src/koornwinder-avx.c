/* This file is part of SQT, a library for Singular Quadrature on Triangles
 *
 * Copyright (C) 2020, 2021 Michael Carley
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <glib.h>

#include <sqt.h>

#include "sqt-private.h"

#ifdef HAVE_AVX_INSTRUCTIONS
#include <immintrin.h>
#endif /*HAVE_AVX_INSTRUCTIONS*/

#if defined(SQT_SINGLE_PRECISION) || !defined(HAVE_AVX_INSTRUCTIONS)
static gint koornwinder_deriv_recursion4(gint N, gint m,
					 SQT_REAL *Pm, SQT_REAL *dPm,
					 SQT_REAL *vpm,
					 SQT_REAL *u, SQT_REAL *v,
					 SQT_REAL *w, SQT_REAL *vinv,
					 SQT_REAL *K, gint kstr, gint offk,
					 SQT_REAL *Ku, gint ustr, gint offu,
					 SQT_REAL *Kv, gint vstr, gint offv,
					 gint nst)

{
  SQT_REAL *tmp, *Jnm1, *Jnm, buf[8], dJnm[4] ;
  gint n, idx, i ;

  Jnm1 = &(buf[0]) ; Jnm = &(Jnm1[4]) ;
  for ( i = 0 ; i < 4 ; i ++ ) {
    /* dPm[i] *= 2.0 ; */
    Jnm1[i] = 1.0 ; dJnm[i] = 0.0 ;
  }
  n = m ;
  idx = n*(n+1)/2 + m ;

  if ( idx < nst ) {
    for ( i = 0 ; i < 4 ; i ++ ) {
      K [kstr*(i*offk+idx)] = sqrt((1+2*m)*(n+1))*vpm[i] ;
      Ku[ustr*(i*offu+idx)] = dPm[i]*K[kstr*(i*offk+idx)]*vinv[i] ;
      K [kstr*(i*offk+idx)] *= Pm[i] ;
      Kv[vstr*(i*offv+idx)] = Ku[ustr*(i*offu+idx)]*u[i]*vinv[i] -
	K[kstr*(i*offk+idx)]*(m*vinv[i]) ;
    }
  }
  if ( m+1 > N ) return 0 ;
  for ( i = 0 ; i < 4 ; i ++ ) {  
    Jnm[i] = w[i]*(1.5 + m) - 0.5 - m ; dJnm[i] = 1.5 + m ;
  }
  n = m + 1 ;
  idx = n*(n+1)/2 + m ;
  if ( idx < nst ) {
    for ( i = 0 ; i < 4 ; i ++ ) {
      K [kstr*(i*offk+idx)] = sqrt((1+2*m)*(n+1))*vpm[i] ;
      Ku[ustr*(i*offu+idx)] = Jnm[i]*dPm[i]*K[kstr*(i*offk+idx)]*vinv[i] ;
      K [kstr*(i*offk+idx)] *= Pm[i] ;
      Kv[vstr*(i*offv+idx)] =
	Ku[ustr*(i*offu+idx)]*u[i]*vinv[i] -
	K[kstr*(i*offk+idx)]*(Jnm[i]*m*vinv[i] + dJnm[i]*2.0) ;
      K[kstr*(i*offk+idx)] *= Jnm[i] ;
    }
  }

  if ( N*(N+1)/2 + m >= nst ) N -- ;

  for ( n = m+2 ; n <= N ; n ++ ) {
    tmp = Jnm ;
    for ( i = 0 ; i < 4 ; i ++ ) {
      Jnm1[i] = (Jnm[i]*(w[i]*(4*n*n-1) - (2*m+1)*(2*m+1))*n -
		 Jnm1[i]*(n-1-m)*(n+m)*(2*n+1))/((n-m)*(n+m+1)*(2*n-1)) ;
    }
    Jnm = Jnm1 ; Jnm1 = tmp ;
    for ( i = 0 ; i < 4 ; i ++ ) {
      dJnm[i]  = (2.0*(n+m+1)*Jnm1[i] - ((2*n+1)*w[i] + 2*m+1)*Jnm[i])*(n-m) ;
      dJnm[i] /= (1.0-w[i]*w[i])*(2*n+1) ;
    }

    idx = n*(n+1)/2 + m ;
    for ( i = 0 ; i < 4 ; i ++ ) {
      K [kstr*(i*offk+idx)] = sqrt((1+2*m)*(n+1))*vpm[i] ;
      Ku[ustr*(i*offu+idx)] = Jnm[i]*dPm[i]*K[kstr*(i*offk+idx)]*vinv[i] ;
      K [kstr*(i*offk+idx)] *= Pm[i] ;
      Kv[vstr*(i*offv+idx)] =
	Ku[ustr*(i*offu+idx)]*u[i]*vinv[i] -
	K[kstr*(i*offk+idx)]*(Jnm[i]*m*vinv[i] + dJnm[i]*2.0) ;
      K[kstr*(i*offk+idx)] *= Jnm[i] ;
    }
  }
  
  return 0 ;
}
#else
static gint koornwinder_deriv_recursion4(gint N, gint m,
					 gdouble *Pm, gdouble *dPm,
					 gdouble *vpm,
					 gdouble *u, gdouble *v,
					 gdouble *w, gdouble *vinv,
					 gdouble *K, gint kstr, gint offk,
					 gdouble *Ku, gint ustr, gint offu,
					 gdouble *Kv, gint vstr, gint offv,
					 gint nst)

{
  gdouble Kb0[4], Kb1[4], Kb2[4], cnm ;
  __m256d rPm, rJnm, rJnm1, rdJnm, op1, op2, op3, rcnm, rw, rvpm, rvinv ;
  __m256d rm, rdPm, ruvinv ;
  gint n, idx, i ;

  rPm   = _mm256_loadu_pd(Pm) ;
  rdPm  = _mm256_loadu_pd(dPm) ;
  /* rdPm  = _mm256_mul_pd(rdPm, _mm256_set1_pd(2.0)) ; */
  rw    = _mm256_loadu_pd(w) ;
  rm    = _mm256_set1_pd((gdouble)m) ;
  rvpm  = _mm256_loadu_pd(vpm) ;
  rvinv = _mm256_loadu_pd(vinv) ;
  ruvinv = _mm256_loadu_pd(u) ; ruvinv = _mm256_mul_pd(ruvinv, rvinv) ;

  rJnm1 = _mm256_set1_pd(1.0) ;
  /* rdJnm = _mm256_set1_pd(0.0) ; */

  n = m ;
  idx = n*(n+1)/2 + m ;
  cnm = (1+2*m)*(n+1) ;
  rcnm = _mm256_set1_pd(cnm) ; rcnm = _mm256_sqrt_pd(rcnm) ;

  op1 = _mm256_mul_pd(rcnm, rvpm) ;
  op2 = _mm256_mul_pd(rdPm, op1) ;
  op2 = _mm256_mul_pd(op2, rvinv) ;
  op1 = _mm256_mul_pd(op1, rPm) ;
  _mm256_storeu_pd(Kb0, op1) ;
  _mm256_storeu_pd(Kb1, op2) ;
  /* op2 = _mm256_mul_pd(op2,ru) ; op2 = _mm256_mul_pd(op2,rvinv) ; */
  op2 = _mm256_mul_pd(op2,ruvinv) ;
  op1 = _mm256_mul_pd(op1, rvinv) ; op1 = _mm256_mul_pd(op1, rm) ;
  op1 = _mm256_sub_pd(op2, op1) ;
  _mm256_storeu_pd(Kb2, op1) ;
  if ( idx < nst ) {
    for ( i = 0 ; i < 4 ; i ++ ) {
      K [kstr*(i*offk+idx)] = Kb0[i] ;
      Ku[ustr*(i*offu+idx)] = Kb1[i] ;
      Kv[vstr*(i*offv+idx)] = Kb2[i] ;
    }
  }
  /* if ( m+1 > N ) return 0 ; */

  cnm = 1.5 + m ; rdJnm = _mm256_set1_pd(cnm) ;
  rJnm = _mm256_mul_pd(rdJnm, rw) ;
  rJnm = _mm256_sub_pd(rJnm, _mm256_set1_pd(m+0.5)) ;
  rdJnm = _mm256_mul_pd(rdJnm, _mm256_set1_pd(2.0)) ;
  
  n = m + 1 ;
  idx = n*(n+1)/2 + m ;
  cnm = (1+2*m)*(n+1) ;
  rcnm = _mm256_set1_pd(cnm) ; rcnm = _mm256_sqrt_pd(rcnm) ;

  op1 = _mm256_mul_pd(rcnm, rvpm) ;    /*cnm*vpm*/
  op2 = _mm256_mul_pd(rdPm, op1) ;
  op2 = _mm256_mul_pd(op2, rvinv) ;
  op2 = _mm256_mul_pd(op2, rJnm) ;
  op1 = _mm256_mul_pd(op1, rPm) ;
  op3 = _mm256_mul_pd(op1, rJnm) ;     /*final K*/
  _mm256_storeu_pd(Kb0, op3) ;
  _mm256_storeu_pd(Kb1, op2) ;
  /* op2 = _mm256_mul_pd(op2,ru) ; op2 = _mm256_mul_pd(op2,rvinv) ; */
  op2 = _mm256_mul_pd(op2,ruvinv) ;
  op3 = _mm256_mul_pd(rJnm, rm) ; op3 = _mm256_mul_pd(op3, rvinv) ;
  op3 = _mm256_add_pd(op3, rdJnm) ;
  op1 = _mm256_mul_pd(op1, op3) ;
  op1 = _mm256_sub_pd(op2, op1) ;
  _mm256_storeu_pd(Kb2, op1) ;
  if ( idx < nst ) {
    for ( i = 0 ; i < 4 ; i ++ ) {
      K [kstr*(i*offk+idx)] = Kb0[i] ;
      Ku[ustr*(i*offu+idx)] = Kb1[i] ;
      Kv[vstr*(i*offv+idx)] = Kb2[i] ;
    }
  }

  if ( N*(N+1)/2 + m >= nst ) N -- ;

  for ( n = m+2 ; n <= N ; n ++ ) {
    op1 = _mm256_mul_pd(rw, _mm256_set1_pd((gdouble)(4*n*n-1)*n)) ;
    op1 = _mm256_sub_pd(op1, _mm256_set1_pd((gdouble)(2*m+1)*(2*m+1)*n)) ;
    op1 = _mm256_mul_pd(op1, rJnm) ;
    op2 = _mm256_mul_pd(rJnm1,
			_mm256_set1_pd((gdouble)((n-1-m)*(n+m)*(2*n+1)))) ;
    op1 = _mm256_sub_pd(op1, op2) ;
    op1 = _mm256_div_pd(op1,
			_mm256_set1_pd((gdouble)((n-m)*(n+m+1)*(2*n-1)))) ;
    rJnm1 = rJnm ; rJnm = op1 ;
    op1 = _mm256_mul_pd(rJnm1, _mm256_set1_pd(2.0*(n+m+1))) ;
    op2 = _mm256_mul_pd(rw,  _mm256_set1_pd(2.0*n+1)) ;
    op2 = _mm256_add_pd(op2, _mm256_set1_pd(2.0*m+1)) ;
    op2 = _mm256_mul_pd(op2, rJnm) ;
    op1 = _mm256_sub_pd(op1, op2) ;
    op1 = _mm256_mul_pd(op1, _mm256_set1_pd(2.0*(n-m))) ;
    op2 = _mm256_mul_pd(rw, rw) ;
    op2 = _mm256_sub_pd(_mm256_set1_pd(1.0), op2) ;
    op2 = _mm256_mul_pd(op2, _mm256_set1_pd(2.0*n+1)) ;
    rdJnm = _mm256_div_pd(op1, op2) ;
    
    idx = n*(n+1)/2 + m ;
    cnm = (1+2*m)*(n+1) ;
    rcnm = _mm256_set1_pd(cnm) ; rcnm = _mm256_sqrt_pd(rcnm) ;

    op1 = _mm256_mul_pd(rcnm, rvpm) ;    /*cnm*vpm*/
    op2 = _mm256_mul_pd(rdPm, op1) ;
    op2 = _mm256_mul_pd(op2, rvinv) ;
    op2 = _mm256_mul_pd(op2, rJnm) ;
    op1 = _mm256_mul_pd(op1, rPm) ;
    op3 = _mm256_mul_pd(op1, rJnm) ;     /*final K*/
    _mm256_storeu_pd(Kb0, op3) ;
    _mm256_storeu_pd(Kb1, op2) ;
    /* op2 = _mm256_mul_pd(op2,ru) ; op2 = _mm256_mul_pd(op2,rvinv) ; */
    op2 = _mm256_mul_pd(op2,ruvinv) ;
    op3 = _mm256_mul_pd(rJnm, rm) ; op3 = _mm256_mul_pd(op3, rvinv) ;
    op3 = _mm256_add_pd(op3, rdJnm) ;
    op1 = _mm256_mul_pd(op1, op3) ;
    op1 = _mm256_sub_pd(op2, op1) ;
    _mm256_storeu_pd(Kb2, op1) ;
    for ( i = 0 ; i < 4 ; i ++ ) {
      K [kstr*(i*offk+idx)] = Kb0[i] ;
      Ku[ustr*(i*offu+idx)] = Kb1[i] ;
      Kv[vstr*(i*offv+idx)] = Kb2[i] ;
    }
  }
  
  return 0 ;
}
#endif /*defined(SQT_SINGLE_PRECISION) || !defined(HAVE_AVX_INSTRUCTIONS)*/

gint SQT_FUNCTION_NAME(sqt_koornwinder_deriv_nm_vector)(gint N,
							SQT_REAL *u,
							SQT_REAL *v,
							gint nu,
							SQT_REAL *K ,
							gint kstr, gint offk,
							SQT_REAL *Ku,
							gint ustr, gint offu,
							SQT_REAL *Kv,
							gint vstr, gint offv,
							gint nst)

/*
  Definition of Koornwinder polynomials from Greengard et al, Fast
  multipole methods for the evaluation of layer potentials with
  locally-corrected quadratures arXiv:2006.02545v1.

  Jacobi and Legendre recursions from DL2MF

  indexing idx_{nm} = n*(n+1)/2 + m
*/
  
{
  gint m, i ;
  SQT_REAL x[4], *tmp, w[4], *Pm, *Pmm1, dPm[4], vpm[4], vinv[4], buf[8] ;

  g_assert(nu <= 4) ;
  for ( i = 0 ; i < 4 ; i ++ ) {
    w[i] = 1.0 - 2.0*v[i] ; vinv[i] = 1.0/(1.0 - v[i]) ;
    vpm[i] = (1.0 - v[i]) ;
    x[i] = (2.0*u[i] + -vpm[i])*vinv[i] ;
    vpm[i] *= M_SQRT2 ;
  }
  
  m = 0 ;
  Pm = &(buf[0]) ;
  for ( i = 0 ; i < 4 ; i ++ ) {
    Pm[i] = 1.0 ; dPm[i] = 0.0 ;
    buf[4+i] = M_SQRT2 ;
  }
  koornwinder_deriv_recursion4(N, m, Pm, dPm, &(buf[4]), u, v, w, vinv,
			       K, kstr, offk, Ku, ustr, offu, Kv, vstr, offv,
			       nst) ;

  m = 1 ;
  Pmm1 = Pm ; Pm = &(Pmm1[4]) ;
  for ( i = 0 ; i < 4 ; i ++ ) {
    /*derivative of Pm scaled by 2.0 to save a multiplication in the
      recursion*/
    Pm[i] = x[i] ; dPm[i] = 2.0 ;
  }

  koornwinder_deriv_recursion4(N, m, Pm, dPm, vpm, u, v, w, vinv,
			       K, kstr, offk, Ku, ustr, offu, Kv, vstr, offv,
			       nst) ;

  /*recursive generation of Legendre polynomials*/
  for ( m = 2 ; m <= N ; m ++ ) {
    tmp = Pm ;
    for ( i = 0 ; i < 4 ; i ++ ) {
      Pmm1[i] = (Pm[i]*x[i]*(2*m-1) - Pmm1[i]*(m-1))/m ;
    }
    Pm = Pmm1 ; Pmm1 = tmp ;
    for ( i = 0 ; i < 4 ; i ++ ) {
      /*derivative of Pm scaled by 2.0 to save a multiplication in the
	recursion*/
      dPm[i] = 2.0*(Pmm1[i] - x[i]*Pm[i])/(1.0-x[i]*x[i])*m ;
      vpm[i] *= 1.0 - v[i] ;
    }
    koornwinder_deriv_recursion4(N, m, Pm, dPm, vpm, u, v, w, vinv,
				 K, kstr, offk, Ku, ustr, offu, Kv, vstr, offv,
				 nst) ;
  }
  
  return 0 ;
}

#if defined(SQT_SINGLE_PRECISION) || !defined(HAVE_AVX_INSTRUCTIONS)
static gint koornwinder_recursion4(gint N, gint m,
				   SQT_REAL *Pm, SQT_REAL *w,
				   SQT_REAL *Knm, gint str, gint nst)

{
  SQT_REAL cnm, *tmp, *Jnm1, *Jnm, buf[8], A, B, C ;
  gint n, idx, i ;

  Jnm1 = &(buf[0]) ; Jnm = &(Jnm1[4]) ;
  for ( i = 0 ; i < 4 ; i ++ ) Jnm1[i] = 1.0 ;
  n = m ;
  idx = n*(n+1)/2 + m ;
  cnm = sqrt(2.0*(1+2*m)*(m+1)) ;

  if ( idx < nst ) {
    for ( i = 0 ; i < 4 ; i ++ ) 
      Knm[(i*nst+idx)*str] = Jnm1[i]*Pm[i]*cnm ;
  }
  if ( m+1 > N ) return 0 ;
  
  A = 0.5*(2*m+3) ; B = -0.5*(2*m+1) ;

  for ( i = 0 ; i < 4 ; i ++ ) Jnm[i] = A*w[i] + B ;
  n = m + 1 ;
  cnm = sqrt(2.0*(1+2*m)*(n+1)) ;
  idx = n*(n+1)/2 + m ;
  if ( idx < nst ) {
    for ( i = 0 ; i < 4 ; i ++ ) 
      Knm[(i*nst+idx)*str] = Jnm[i]*Pm[i]*cnm ;
  }

  if ( N*(N+1)/2 + m >= nst ) N -- ;

  for ( n = m+2 ; n <= N ; n ++ ) {
    A =  (SQT_REAL)(n-1+1)*(2*n-2+3)/(n-1-m+1)/(n-1+m+2) ;
    B = -(SQT_REAL)(2*m+1)*(2*m+1)*(n-1+1)/(n-1-m+1)/(n-1+m+2)/(2*n-2+1) ;
    C =  (SQT_REAL)(n-1-m)*(n-1+m+1)*(2*n-2+3)/(n-1-m+1)/(n-1+m+2)/(2*n-2+1) ;
    
    tmp = Jnm ;
    for ( i = 0 ; i < 4 ; i ++ ) {
      Jnm1[i] = (A*w[i] + B)*Jnm[i] - C*Jnm1[i] ;
    }
    Jnm = Jnm1 ; Jnm1 = tmp ;
    cnm = sqrt(2.0*(1+2*m)*(n+1)) ;
    idx = n*(n+1)/2 + m ;
    for ( i = 0 ; i < 4 ; i ++ ) Knm[(i*nst+idx)*str] = Jnm[i]*Pm[i]*cnm ;
  }
  
  return 0 ;
}
#else
static gint koornwinder_recursion4(gint N, gint m,
				   gdouble *Pm, gdouble *w,
				   gdouble *Knm, gint str, gint nst)

{
  gdouble cnm, A, B, C, Kb[4] ;
  __m256d rPm, rJnm, rJnm1, op1, op2, rcnm, rw ;
  gint n, idx ;

  rPm = _mm256_loadu_pd(Pm) ;
  rw  = _mm256_loadu_pd(w) ;

  rJnm1 = _mm256_set1_pd(1.0) ;
  n = m ;
  idx = n*(n+1)/2 + m ;
  cnm = 2.0*(1+2*m)*(m+1) ;
  rcnm = _mm256_set1_pd(cnm) ; rcnm = _mm256_sqrt_pd(rcnm) ;

  op1 = _mm256_mul_pd(rJnm1,rPm) ;
  op1 = _mm256_mul_pd(op1,rcnm) ;
  _mm256_storeu_pd(Kb, op1) ;
  if ( idx < nst ) {
    Knm[(0*nst+idx)*str] = Kb[0] ;
    Knm[(1*nst+idx)*str] = Kb[1] ;
    Knm[(2*nst+idx)*str] = Kb[2] ;
    Knm[(3*nst+idx)*str] = Kb[3] ;
  }
  if ( m+1 > N ) return 0 ;
  
  A = 0.5*(2*m+3) ; B = -0.5*(2*m+1) ;
  op1 = _mm256_set1_pd(A) ; op2 = _mm256_set1_pd(B) ;
  rJnm = _mm256_mul_pd(rw, op1) ;
  rJnm = _mm256_add_pd(rJnm, op2) ;
  
  n = m + 1 ;
  cnm = 2.0*(1+2*m)*(n+1) ;
  rcnm = _mm256_set1_pd(cnm) ; rcnm = _mm256_sqrt_pd(rcnm) ;
  idx = n*(n+1)/2 + m ;
  op1 = _mm256_mul_pd(rJnm,rPm) ;
  op1 = _mm256_mul_pd(op1,rcnm) ;
  _mm256_storeu_pd(Kb, op1) ;

  if ( idx < nst ) {
    Knm[(0*nst+idx)*str] = Kb[0] ;
    Knm[(1*nst+idx)*str] = Kb[1] ;
    Knm[(2*nst+idx)*str] = Kb[2] ;
    Knm[(3*nst+idx)*str] = Kb[3] ;
  }

  if ( N*(N+1)/2 + m >= nst ) N -- ;

  for ( n = m+2 ; n <= N ; n ++ ) {
    cnm = 2.0*(1+2*m)*(n+1) ;
    rcnm = _mm256_set1_pd(cnm) ; rcnm = _mm256_sqrt_pd(rcnm) ;
    A =  (gdouble)(n-1+1)*(2*n-2+3)/(n-1-m+1)/(n-1+m+2) ;
    B = -(gdouble)(2*m+1)*(2*m+1)*(n-1+1)/(n-1-m+1)/(n-1+m+2)/(2*n-2+1) ;
    C =  (gdouble)(n-1-m)*(n-1+m+1)*(2*n-2+3)/(n-1-m+1)/(n-1+m+2)/(2*n-2+1) ;
    
    op1 = _mm256_set1_pd(A) ;
    op2 = _mm256_set1_pd(B) ;
    op1 = _mm256_mul_pd(rw, op1) ;      /*A*w*/
    op1 = _mm256_add_pd(op1, op2) ;     /*A*w + B*/
    op1 = _mm256_mul_pd(op1, rJnm) ;    /*(A*w + B)*Jnm*/
    op2 = _mm256_set1_pd(C) ;
    op2 = _mm256_mul_pd(op2, rJnm1) ;   /*C*Jnm1*/
    op1 = _mm256_sub_pd(op1, op2) ;     /*(A*w + B)*Jnm - C*Jnm1*/
    rJnm1 = rJnm ; rJnm = op1 ;
    op1 = _mm256_mul_pd(rJnm,rPm) ;
    op1 = _mm256_mul_pd(op1,rcnm) ;
    _mm256_storeu_pd(Kb, op1) ;

    idx = n*(n+1)/2 + m ;
    Knm[(0*nst+idx)*str] = Kb[0] ;
    Knm[(1*nst+idx)*str] = Kb[1] ;
    Knm[(2*nst+idx)*str] = Kb[2] ;
    Knm[(3*nst+idx)*str] = Kb[3] ;
  }

  return 0 ;
}
#endif /*defined(SQT_SINGLE_PRECISION) || !defined(HAVE_AVX_INSTRUCTIONS)*/

gint SQT_FUNCTION_NAME(sqt_koornwinder_nm_vector)(gint N,
						  SQT_REAL *u, SQT_REAL *v,
						  gint nu,
						  SQT_REAL *Knm, gint str,
						  gint nst)

/*
  Definition of Koornwinder polynomials from Greengard et al, Fast
  multipole methods for the evaluation of layer potentials with
  locally-corrected quadratures arXiv:2006.02545v1.

  Jacobi and Legendre recursions from DLMF

  indexing idx_{nm} = n*(n+1)/2 + m
*/
  
{
  gint m, i ;
  SQT_REAL x[4], *tmp, w[4], *Pm, *Pmm1, buf[8] ;

  g_assert(nu <= 4) ;

  Pm = &(buf[0]) ; 
  for ( i = 0 ; i < 4 ; i ++ ) {
    w[i] = 1.0 - 2.0*v[i] ;
    x[i] = (2.0*u[i] + v[i] - 1.0)/(1.0 - v[i]) ;
    Pm[i] = 1.0 ;
    Pm[4+i] = (1.0 - v[i])*x[i] ;
  }
  m = 0 ;

  koornwinder_recursion4(N, m, Pm, w, Knm, str, nst) ;

  m = 1 ;
  Pmm1 = Pm ; Pm = &(buf[4]) ;
  koornwinder_recursion4(N, m, Pm, w, Knm, str, nst) ;

  /*recursive generation of Legendre polynomials and Jacobi recursion*/
  for ( m = 2 ; m <= N ; m ++ ) {
    tmp = Pm ;
    for ( i = 0 ; i < 4 ; i ++ ) {
      Pmm1[i] = (1.0-v[i])*((2.0*m-1)*x[i]*Pm[i] -
			    (1.0-v[i])*(m-1)*Pmm1[i])/m ;
    }
    Pm = Pmm1 ; Pmm1 = tmp ;
    koornwinder_recursion4(N, m, Pm, w, Knm, str, nst) ;
  }  
  
  return 0 ;
}
