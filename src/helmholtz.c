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

#include "config.h"

#include "sqt-private.h"

#ifdef HAVE_AVX_INSTRUCTIONS
#include <immintrin.h>
#endif /*HAVE_AVX_INSTRUCTIONS*/

#if 0

static gint laplace_quad_weights_tri(SQT_REAL s, SQT_REAL t,
				     SQT_REAL w,
				     SQT_REAL *y, SQT_REAL *n,
				     SQT_REAL *K, gint nk,
				     SQT_REAL *quad, gint nc,
				     gint init,
				     gpointer data[])

{
  SQT_REAL *Knm = data[SQT_DATA_KNM] ;
  gint nq = *((gint *)(data[SQT_DATA_NKNM])) ;
  gint Nk = *((gint *)(data[SQT_DATA_ORDER_K])) ;
  SQT_REAL *x  = data[SQT_DATA_TARGET] ;
  SQT_REAL *Kq  = data[SQT_DATA_MATRIX] ;
  SQT_REAL R, G, dG, wt, d1 = 1.0 ;
  gint i1 = 1 ;

  /*Koornwinder polynomials at evaluation point*/
  SQT_FUNCTION_NAME(sqt_koornwinder_nm)(Nk, s, t, Knm, 1, nq) ;

  R = sqt_vector_distance(x, y) ;

  G = 0.25*M_1_PI/R ;

  dG = sqt_vector_diff_scalar(x, y, n)/R/R*G ;

#ifndef SQT_SINGLE_PRECISION  
  wt = w*G ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   0]), i1) ;
  wt = w*dG ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[nc/2]), i1) ;
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ;
  wt = w*G ;
  blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   0]), i1) ;
  wt = w*dG ;
  blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[nc/2]), i1) ;
#endif
  
  return 0 ;
}

static gint laplace_quad_weights_cached(SQT_REAL s, SQT_REAL t,
					SQT_REAL w,
					SQT_REAL *y, SQT_REAL *n,
					SQT_REAL *K, gint nk,
					SQT_REAL *quad, gint nc,
					gint init,
					gpointer data[])

{
  SQT_REAL *Knm = &(y[SQT_CACHE_STRIDE]) ;
  gint nq = *((gint *)(data[SQT_DATA_NKNM])) ;
  gint Nk = *((gint *)(data[SQT_DATA_ORDER_K])) ;
  SQT_REAL *x  = data[SQT_DATA_TARGET] ;
  SQT_REAL *Kq  = data[SQT_DATA_MATRIX] ;
  SQT_REAL R, G, dG, wt, d1 = 1.0 ;
  gint i1 = 1 ;

  /*Koornwinder polynomials at evaluation point*/
  if ( init != 0 )
    SQT_FUNCTION_NAME(sqt_koornwinder_nm)(Nk, s, t, Knm, 1, nq) ;
  R = sqt_vector_distance(x, y) ;

  G = 0.25*M_1_PI/R ;

  dG = sqt_vector_diff_scalar(x, y, n)/R/R*G ;

#ifndef SQT_SINGLE_PRECISION  
  wt = w*G ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   0]), i1) ;
  wt = w*dG ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[nc/2]), i1) ;
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ;
  wt = w*G ;
  blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   0]), i1) ;
  wt = w*dG ;
  blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[nc/2]), i1) ;
#endif
  
  return 0 ;
}

static gint laplace_quad_matrix(SQT_REAL s, SQT_REAL t,
				SQT_REAL w,
				SQT_REAL *y, SQT_REAL *n,
				SQT_REAL *Knm, gint nk,
				SQT_REAL *quad, gint nc,
				gint init,
				gpointer data[])

{
  gint nq = *((gint *)(data[SQT_DATA_NKNM])) ;
  SQT_REAL *x  = data[SQT_DATA_TARGET] ;
  gint xstr = *((gint *)data[SQT_DATA_STRIDE]) ;
  gint nx = *((gint *)data[SQT_DATA_NUMBER]) ;
  SQT_REAL *Kq  = data[SQT_DATA_MATRIX] ;
  SQT_REAL G, dG, d0 = 0.0 ;
  SQT_REAL work[453], r[3] ;
  gint i1 = 1, i, ns ;

  G = w*0.25*M_1_PI ;
#ifndef SQT_SINGLE_PRECISION  
  blaswrap_dgemv(TRUE, nq, nq, G, Kq, nq, Knm, i1, d0, work, i1) ;
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ;
  blaswrap_sgemv(TRUE, nq, nq, G, Kq, nq, Knm, i1, d0, work, i1) ;
#endif /*SQT_SINGLE_PRECISION*/
  
  ns = nc/nx ;
  for ( i = 0 ; i < nx ; i ++ ) {
    sqt_vector_diff(r, &(x[i*xstr]), y) ;
    dG = sqt_vector_scalar(r, n) ;
    G = 1.0/sqt_vector_length(r) ;
    dG *= G*G*G ;
    
#ifndef SQT_SINGLE_PRECISION  
    blaswrap_daxpy(nq,  G, work, i1, &(quad[i*ns+   0]), i1) ;
    blaswrap_daxpy(nq, dG, work, i1, &(quad[i*ns+ns/2]), i1) ;
#else /*SQT_SINGLE_PRECISION*/
    g_assert_not_reached() ;
    blaswrap_saxpy(nq,  G, work, i1, &(quad[i*ns+   0]), i1) ;
    blaswrap_saxpy(nq, dG, work, i1, &(quad[i*ns+ns/2]), i1) ;
#endif /*SQT_SINGLE_PRECISION*/
  }  

  /*this works but doesn't give a very large speed up*/
/* #ifndef SQT_SINGLE_PRECISION */
/*   gdouble al, bt ; */
/*   al = 1.0 ; bt = 1.0 ; */
/*   blaswrap_dgemm(FALSE, FALSE, nx, nq, i1, al,  g, i1, work, nq, bt, */
/* 		 &(quad[0*ns/2]), ns) ; */
/*   blaswrap_dgemm(FALSE, FALSE, nx, nq, i1, al, dg, i1, work, nq, bt, */
/* 		 &(quad[1*ns/2]), ns) ; */
/* #else /\*SQT_SINGLE_PRECISION*\/ */
/*     /\*this is to keep the compiler quiet about unused variables*\/ */
/*     /\* wt = 0.0 ; Kq[0] = sin(wt*d1*i1*dG) ; *\/ */
/*     g_assert_not_reached() ; */
/* #endif /\*SQT_SINGLE_PRECISION*\/ */

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_laplace_weights_tri_basic)(SQT_REAL *xe,
						      gint xstr, gint ne,
						      SQT_REAL *q, gint nq,
						      SQT_REAL *Kq,
						      gint nqk,
						      gint nK,
						      SQT_REAL *x,
						      SQT_REAL *w)

{
  SQT_REAL Knm[2048] ;
  gpointer data[SQT_DATA_WIDTH] ;
  
  data[SQT_DATA_ELEMENT]   = xe ; 
  data[SQT_DATA_STRIDE]    = &xstr ;
  data[SQT_DATA_NUMBER]    = &(ne) ;
  data[SQT_DATA_TARGET]    = x ; 
  data[SQT_DATA_MATRIX]    = Kq ;
  data[SQT_DATA_KNM]       = Knm ;
  data[SQT_DATA_NKNM]      = &nqk ;
  data[SQT_DATA_ORDER_K]   = &nK ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)laplace_quad_weights_tri ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights_tri ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_basic_quad_tri)(xe, xstr, ne, q, nq, func,
					w, 2*nqk, data) ;
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_laplace_weights_tri_adaptive)(SQT_REAL *xe,
							 gint xstr, gint ne,
							 SQT_REAL *q, gint nq,
							 SQT_REAL *Kq,
							 gint nqk,
							 gint nK,
							 SQT_REAL tol,
							 gint dmax,
							 SQT_REAL *x,
							 SQT_REAL *w)

{
  SQT_REAL Knm[2048] ;
  gpointer data[SQT_DATA_WIDTH] ;
  
  data[SQT_DATA_ELEMENT]   = xe ; 
  data[SQT_DATA_STRIDE]    = &xstr ;
  data[SQT_DATA_NUMBER]    = &(ne) ;
  data[SQT_DATA_TARGET]    = x ; 
  data[SQT_DATA_MATRIX]    = Kq ;
  data[SQT_DATA_KNM]       = Knm ;
  data[SQT_DATA_NKNM]      = &nqk ;
  data[SQT_DATA_ORDER_K]   = &nK ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)laplace_quad_weights_tri ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights_tri ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_adaptive_quad_tri)(xe, xstr, ne, q, nq, func,
					   w, 2*nqk, tol, dmax, data) ;
  
  return 0 ;
}


gint SQT_FUNCTION_NAME(sqt_laplace_matrix_kw_adaptive)(SQT_REAL *ce,
						       gint ne, gint Nk,
						       SQT_REAL *Kq,
						       SQT_REAL *q, gint nq,
						       SQT_REAL tol,
						       gint dmax,
						       SQT_REAL *x,
						       gint xstr, gint nx,
						       SQT_REAL *Ast,
						       SQT_REAL *work)

{
  gpointer data[SQT_DATA_WIDTH] ;
  
  data[SQT_DATA_TARGET]    = x ; 
  data[SQT_DATA_STRIDE]    = &xstr ; 
  data[SQT_DATA_NUMBER]    = &nx ;
  data[SQT_DATA_MATRIX]    = Kq ;
  data[SQT_DATA_NKNM]      = &ne ;
  data[SQT_DATA_ORDER_K]   = &Nk ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)laplace_quad_matrix ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_matrix ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_adaptive_quad_kw)(ce, ne, Nk, q, nq, func,
					  Ast, 2*ne*nx, tol, dmax,
					  data, work) ;
  
  return 0 ;
}


gint SQT_FUNCTION_NAME(sqt_laplace_weights_tri_singular)(SQT_REAL *xe,
							 gint xstr, gint ne,
							 SQT_REAL *Kq,
							 gint nqk,
							 gint nK,
							 gint N,
							 SQT_REAL s0,
							 SQT_REAL t0,
							 SQT_REAL *w)

{
  SQT_REAL Knm[512], x[3], n[3], J ;
  gpointer data[SQT_DATA_WIDTH] ;
  
  SQT_FUNCTION_NAME(sqt_element_point_3d)(xe, xstr, ne, s0, t0, x, n, &J) ;

  data[SQT_DATA_ELEMENT]   = xe ; 
  data[SQT_DATA_STRIDE]    = &xstr ;
  data[SQT_DATA_NUMBER]    = &(ne) ;
  data[SQT_DATA_TARGET]    = x ; 
  data[SQT_DATA_MATRIX]    = Kq ;
  data[SQT_DATA_KNM]       = Knm ;
  data[SQT_DATA_NKNM]      = &nqk ;
  data[SQT_DATA_ORDER_K]   = &nK ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)laplace_quad_weights_tri ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights_tri ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_singular_quad_tri)(xe, xstr, ne, s0, t0, N, func,
					   w, 2*nqk, data) ;
  
  return 0 ;
}

gint
SQT_FUNCTION_NAME(sqt_laplace_source_target_tri_basic)(SQT_REAL *xse,
						       gint xsstr, gint nse,
						       SQT_REAL *q, gint nq,
						       SQT_REAL *Kq,
						       gint nqk, gint nK,
						       SQT_REAL *xte,
						       gint xtstr, gint nte,
						       SQT_REAL *s,
						       gint sstr,
						       SQT_REAL *t,
						       gint tstr,
						       gint nt,
						       SQT_REAL *Ast)
{
  gint i ;
  SQT_REAL x[3], n[3], J ;

  memset(Ast, 0, 2*nqk*nt*sizeof(SQT_REAL)) ;
  for ( i = 0 ; i < nt ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_element_point_3d)(xte, xtstr, nte,
					    s[i*sstr], t[i*tstr], x, n, &J) ;
    SQT_FUNCTION_NAME(sqt_laplace_weights_tri_basic)(xse, xsstr, nse,
						     q, nq, Kq, nqk, nK,
						     x,	&(Ast[i*2*nqk])) ;
  }
  
  return 0 ;
}

gint
SQT_FUNCTION_NAME(sqt_laplace_source_target_tri_adaptive)(SQT_REAL *xse,
							  gint xsstr, gint nse,
							  SQT_REAL *q, gint nq,
							  SQT_REAL *Kq,
							  gint nqk, gint nK,
							  SQT_REAL tol,
							  gint dmax,
							  SQT_REAL *xte,
							  gint xtstr, gint nte,
							  SQT_REAL *s,
							  gint sstr,
							  SQT_REAL *t,
							  gint tstr,
							  gint nt,
							  SQT_REAL *Ast)
{
  gint i ;
  SQT_REAL x[3], n[3], J ;

  memset(Ast, 0, 2*nqk*nt*sizeof(SQT_REAL)) ;
  for ( i = 0 ; i < nt ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_element_point_3d)(xte, xtstr, nte,
					    s[i*sstr], t[i*tstr], x, n, &J) ;
    SQT_FUNCTION_NAME(sqt_laplace_weights_tri_adaptive)(xse, xsstr, nse,
							q, nq, Kq, nqk, nK,
							tol, dmax, x,
							&(Ast[i*2*nqk])) ;
  }
  
  return 0 ;
}

gint
SQT_FUNCTION_NAME(sqt_laplace_source_target_kw_adaptive)(SQT_REAL *xse,
							 gint sstr, gint nse,
							 SQT_REAL *q, gint nq,
							 SQT_REAL *Ks,
							 gint Ns,
							 SQT_REAL tol,
							 gint dmax,
							 SQT_REAL *xte,
							 gint tstr, gint nte,
							 SQT_REAL *Ast,
							 SQT_REAL *work)

/*
 * work space size: 4*dmax*2*nse*nte
 */

{
  gint i3 = 3 ;
  SQT_REAL ce[3*453], al, bt ;

  al = 1.0 ; bt = 0.0 ;
#ifndef SQT_SINGLE_PRECISION
  /*element geometric interpolation coefficients*/
  blaswrap_dgemm(FALSE, FALSE, nse, i3, nse, al, Ks, nse,
		 xse, sstr, bt, ce, i3) ;
#else
  g_assert_not_reached() ; 
  blaswrap_sgemm(FALSE, FALSE, nse, i3, nse, al, Ks, nse,
		 xse, sstr, bt, ce, i3) ;
#endif /*SQT_SINGLE_PRECISION*/
  
  memset(Ast, 0, 2*nse*nte*sizeof(SQT_REAL)) ;  
  SQT_FUNCTION_NAME(sqt_laplace_matrix_kw_adaptive)(ce, nse, Ns, Ks,
  						    q, nq,
  						    tol, dmax,
  						    xte, tstr, nte,
  						    Ast, work) ;  
  return 0 ;
}

static gint laplace_quad_matrix_indexed(SQT_REAL s, SQT_REAL t,
					SQT_REAL w,
					SQT_REAL *y, SQT_REAL *n,
					SQT_REAL *K, gint nk,
					SQT_REAL *quad, gint nc,
					gint init,
					gpointer data[])

{
  SQT_REAL *work = data[SQT_DATA_KNM] ;
  SQT_REAL *x   = data[SQT_DATA_TARGET] ;
  SQT_REAL *Kq  = data[SQT_DATA_MATRIX] ;
  gint nse  = *((gint *)data[SQT_DATA_NKNM]) ;
  /* gint Nk   = *((gint *)data[SQT_DATA_ORDER_K]) ; */
  gint xstr = *((gint *)data[SQT_DATA_STRIDE]) ;
  gint nx   = *((gint *)data[SQT_DATA_NUMBER]) ;
  gint *idx =           data[SQT_DATA_INDICES] ;
  /* gint nq   = *((gint *)data[SQT_DATA_N_QUAD_POINTS]) ; */
  SQT_REAL G, d1 = 1.0, d0 = 0.0, r[3], *g ;
  gint i1 = 1, i, j, ns ;
  gint ncc, na ;
  
  /*Koornwinder polynomials at evaluation point are in input K*/
  G = w*0.25*M_1_PI ;
#ifndef SQT_SINGLE_PRECISION  
  blaswrap_dgemv(TRUE, nse, nse, G, Kq, nse, K, i1, d0, work, i1) ;
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ; /*untested code*/
  blaswrap_sgemv(TRUE, nse, nse, G, Kq, nse, K, i1, d0, work, i1) ;
#endif /*SQT_SINGLE_PRECISION*/

  ns = nc/nx ;
  g = &(work[12*nse]) ;
/* #if defined(SQT_SINGLE_PRECISION) || !defined(HAVE_AVX_INSTRUCTIONS) */
  /*for now this is about as fast as the AVX code below (need to find
    out why)*/
  for ( j = 0 ; j < nx ; j ++ ) {
    SQT_REAL r2 ;
    i = idx[j] ;
    sqt_vector_diff(r, &(x[i*xstr]), y) ;
    g[2*j+1] = sqt_vector_scalar(r, n) ;
    r2 = 1.0/sqt_vector_length2(r) ;
    g[2*j+0] = SQRT(r2) ;
    g[2*j+1] *= r2*g[2*j+0] ;
  }  
  
/* #else /\*defined(SQT_SINGLE_PRECISION) || !defined(HAVE_AVX_INSTRUCTIONS)*\/ */
/*   gint nst  ; */
/*   __m256d rrn, rr, r1, op ; */
/*   r1 = _mm256_set1_pd(1.0) ; */

/*   nst = 4*(nx/4) ; */
/*   for ( j = 0 ; j < nst ; j += 4 ) { */
/*     i = idx[j+0] ; */
/*     sqt_vector_diff(r, &(x[i*xstr]), y) ; */
/*     g[2*j+0] = sqt_vector_length2(r) ; */
/*     g[2*j+4] = sqt_vector_scalar(r, n) ; */
/*     i = idx[j+1] ; */
/*     sqt_vector_diff(r, &(x[i*xstr]), y) ; */
/*     g[2*j+2] = sqt_vector_length2(r) ; */
/*     g[2*j+6] = sqt_vector_scalar(r, n) ; */
/*     i = idx[j+2] ; */
/*     sqt_vector_diff(r, &(x[i*xstr]), y) ; */
/*     g[2*j+1] = sqt_vector_length2(r) ; */
/*     g[2*j+5] = sqt_vector_scalar(r, n) ; */
/*     i = idx[j+3] ; */
/*     sqt_vector_diff(r, &(x[i*xstr]), y) ; */
/*     g[2*j+3] = sqt_vector_length2(r) ; */
/*     g[2*j+7] = sqt_vector_scalar(r, n) ; */
/*     rr  = _mm256_loadu_pd(&(g[2*j+0])) ; /\*R^2*\/ */
/*     rr  = _mm256_div_pd(r1, rr) ; /\*1/R^2*\/ */
/*     rrn = _mm256_loadu_pd(&(g[2*j+4])) ; /\*r.n*\/ */
/*     rrn = _mm256_mul_pd(rrn, rr) ;       /\*r.n/R^2*\/ */
/*     rr  = _mm256_sqrt_pd(rr) ;           /\*1/R*\/ */
/*     rrn = _mm256_mul_pd(rrn, rr) ;       /\*r.n/R^2*\/ */
/*     op  = _mm256_unpacklo_pd(rr, rrn) ; */
/*     _mm256_storeu_pd(&(g[2*j+0]), op) ; */
/*     op  = _mm256_unpackhi_pd(rr, rrn) ; */
/*     _mm256_storeu_pd(&(g[2*j+4]), op) ; */
/*   }   */

/*   for ( j = nst ; j < nx ; j ++ ) { */
/*     gdouble r2 ; */
/*     i = idx[j] ; */
/*     sqt_vector_diff(r, &(x[i*xstr]), y) ; */
/*     g[2*j+1] = sqt_vector_scalar(r, n) ; */
/*     r2 = 1.0/sqt_vector_length2(r) ; */
/*     g[2*j+0] = sqrt(r2) ; */
/*     g[2*j+1] *= r2*g[2*j+0] ; */
/*   } */
  
/* #endif /\*defined(SQT_SINGLE_PRECISION) || !defined(HAVE_AVX_INSTRUCTIONS)*\/ */

  ncc = ns/2 ; na = 2*nx ;
  
#ifndef SQT_SINGLE_PRECISION
  blaswrap_dgemm(FALSE, FALSE, na, nse, i1, d1, g, i1,
		 work, nse, d1, quad, ncc) ;
#else /*SQT_SINGLE_PRECISION  */
  blaswrap_sgemm(FALSE, FALSE, na, nse, i1, d1, g, i1,
		 work, nse, d1, quad, ncc) ;
#endif /*SQT_SINGLE_PRECISION  */

  return 0 ;
}

gint
SQT_FUNCTION_NAME(sqt_laplace_source_indexed_kw_adaptive)(SQT_REAL *xse,
							  gint sstr, gint nse,
							  SQT_REAL *q, gint nq,
							  SQT_REAL *Ks,
							  gint Ns,
							  SQT_REAL tol,
							  gint dmax,
							  SQT_REAL *xte,
							  gint tstr, gint *idx,
							  gint nte,
							  SQT_REAL *Ast,
							  SQT_REAL *work)

/*
 * work space size: 4*dmax*2*nse*nte + 12*nse + 3*nse + 2*nte
 */

{
  gint i3 = 3 ;
  SQT_REAL *ce, *Knm, *awork, al, bt ;
  gpointer data[SQT_DATA_WIDTH] ;

  ce = work ; Knm = &(ce[3*nse]) ; awork = &(Knm[12*nse+2*nte]) ;
  
  al = 1.0 ; bt = 0.0 ;
#ifndef SQT_SINGLE_PRECISION
  /*element geometric interpolation coefficients*/
  blaswrap_dgemm(FALSE, FALSE, nse, i3, nse, al, Ks, nse,
		 xse, sstr, bt, ce, i3) ;
#else
  g_assert_not_reached() ; 
  blaswrap_sgemm(FALSE, FALSE, nse, i3, nse, al, Ks, nse,
		 xse, sstr, bt, ce, i3) ;
#endif /*SQT_SINGLE_PRECISION*/

  data[SQT_DATA_TARGET]    = xte ; 
  data[SQT_DATA_STRIDE]    = &tstr ; 
  data[SQT_DATA_NUMBER]    = &nte ;
  data[SQT_DATA_MATRIX]    = Ks ;
  data[SQT_DATA_KNM]       = Knm ;
  data[SQT_DATA_NKNM]      = &nse ;
  data[SQT_DATA_ORDER_K]   = &Ns ;
  data[SQT_DATA_INDICES]   = idx ;
  data[SQT_DATA_N_QUAD_POINTS] = &nq ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)laplace_quad_matrix_indexed ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_matrix_indexed ;
#endif /*SQT_SINGLE_PRECISION*/

  memset(Ast, 0, 2*nse*nte*sizeof(SQT_REAL)) ;
  SQT_FUNCTION_NAME(sqt_adaptive_quad_kw)(ce, nse, Ns, q, nq, func,
  					  Ast, 2*nse*nte, tol, dmax,
  					  data, awork) ;

  return 0 ;
}

gint
SQT_FUNCTION_NAME(sqt_laplace_source_indexed_kw_cached)(SQT_REAL *xse,
							gint sstr, gint nse,
							SQT_REAL *q, gint nq,
							SQT_REAL *Ks,
							gint Ns,
							SQT_REAL tol,
							gint dmax,
							SQT_REAL *xte,
							gint tstr, gint *idx,
							gint nte,
							gint *icache,
							SQT_REAL *xcache,
							gint cstr,
							SQT_REAL *Ast,
							SQT_REAL *work)

/*
 * work space size: 4*dmax*2*nse*nte
 */

{
  gint i, i3 = 3 ;
  SQT_REAL ce[3*453], al, bt ;
  gpointer data[SQT_DATA_WIDTH] ;
  gboolean init ;

  g_assert(cstr >= SQT_CACHE_STRIDE + nq) ;
  
#ifndef SQT_SINGLE_PRECISION
  /*element geometric interpolation coefficients*/
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, nse, i3, nse, al, Ks, nse,
		 xse, sstr, bt, ce, i3) ;
#else
  g_assert_not_reached() ; 
  al = 1.0 ; bt = 0.0 ;
  blaswrap_sgemm(FALSE, FALSE, nse, i3, nse, al, Ks, nse,
		 xse, sstr, bt, ce, i3) ;
#endif /*SQT_SINGLE_PRECISION*/

  data[SQT_DATA_STRIDE]    = &tstr ; 
  data[SQT_DATA_NUMBER]    = NULL ;
  data[SQT_DATA_MATRIX]    = Ks ;
  data[SQT_DATA_NKNM]      = &nse ;
  data[SQT_DATA_ORDER_K]   = &Ns ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)laplace_quad_weights_cached ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights_cached ;
#endif /*SQT_SINGLE_PRECISION*/

  memset(Ast, 0, 2*nse*nte*sizeof(SQT_REAL)) ;
  init = TRUE ;
  for ( i = 0 ; i < nte ; i ++ ) {
    data[SQT_DATA_TARGET] = &(xte[idx[i]*tstr]) ;
    SQT_FUNCTION_NAME(sqt_cached_quad_kw)(ce, nse, Ns, q, nq, func,
  					  &(Ast[2*i*nse]), 2*nse, tol, dmax,
					  icache, xcache, cstr, init,
					  data, work) ;
    init = FALSE ;
  }

  return 0 ;
}

gint
SQT_FUNCTION_NAME(sqt_laplace_source_target_tri_self)(SQT_REAL *xe,
						      gint xstr, gint ne,
						      SQT_REAL *Kq,
						      gint nqk, gint nK,
						      gint N,
						      SQT_REAL *s,
						      gint sstr,
						      SQT_REAL *t,
						      gint tstr,
						      gint nt,
						      SQT_REAL *Ast)

{
  gint i ;

  memset(Ast, 0, 2*nqk*nqk*sizeof(SQT_REAL)) ;
  for ( i = 0 ; i < nt ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_laplace_weights_tri_singular)(xe, xstr, ne,
							Kq, nqk, nK, N,
							s[i*sstr],
							t[i*tstr],
							&(Ast[i*2*nqk])) ;
  }
  
  return 0 ;
}

#endif

static gint helmholtz_quad_weights(SQT_REAL s, SQT_REAL t,
				   SQT_REAL w,
				   SQT_REAL *y, SQT_REAL *n,
				   SQT_REAL *Knm, gint nk,
				   SQT_REAL *quad, gint nc,
				   gint init,
				   gpointer data[])

{
  gint nq = *((gint *)(data[SQT_DATA_NKNM])) ;
  SQT_REAL *x  = data[SQT_DATA_TARGET] ;
  SQT_REAL *Kq  = data[SQT_DATA_MATRIX] ;
  SQT_REAL k = *((SQT_REAL *)data[SQT_DATA_WAVENUMBER]) ;
  SQT_REAL R, dR, G[2], dG[2], E[2], wt, d1 = 1.0 ;
  gint i1 = 1, i2 = 2 ;

  /*Koornwinder polynomials at evaluation point are in Knm*/
  R = sqt_vector_distance(x, y) ;
  dR = sqt_vector_diff_scalar(x, y, n)/R ;
  E[0] = cos(k*R) ; E[1] = sin(k*R) ;
  
  G[0] = 0.25*M_1_PI/R ;
  G[1] = G[0]*E[1] ;

  dG[0] = (E[0] + k*R*E[1])*G[0]*dR/R ;
  dG[1] = (E[1] - k*R*E[0])*G[0]*dR/R ;

  G[0] *= E[0] ;

#ifndef SQT_SINGLE_PRECISION
  /*it should be possible to do this with one matrix multiplication*/
  wt = w*G[0] ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   0]), i2) ;
  wt = w*G[1] ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   1]), i2) ;
  /*and this*/
  wt = w*dG[0] ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1,
		 &(quad[   nc/2+0]), i2) ;
  wt = w*dG[1] ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1,
		 &(quad[   nc/2+1]), i2) ;
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ; /*untested code*/
  /* wt = w*G ; */
  /* blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   0]), i1) ; */
  /* wt = w*dG ; */
  /* blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[nc/2]), i1) ; */
#endif
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_helmholtz_weights_kw_adaptive)(SQT_REAL *ce,
							  gint ne, gint Nk,
							  SQT_REAL *Kq,
							  SQT_REAL *q, gint nq,
							  SQT_REAL tol,
							  gint dmax,
							  SQT_REAL k,
							  SQT_REAL *x,
							  SQT_REAL *w,
							  SQT_REAL *work)

{
  SQT_REAL Knm[2048] ;
  gpointer data[SQT_DATA_WIDTH] ;
  
  data[SQT_DATA_TARGET]     = x ; 
  data[SQT_DATA_WAVENUMBER] = &k ; 
  data[SQT_DATA_MATRIX]     = Kq ;
  data[SQT_DATA_KNM]        = Knm ;
  data[SQT_DATA_NKNM]       = &ne ;
  data[SQT_DATA_ORDER_K]    = &Nk ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)helmholtz_quad_weights ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)helmholtz_quad_weights ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_adaptive_quad_kw)(ce, ne, Nk, q, nq, func,
					  w, 4*ne, tol, dmax, data, work) ;
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_helmholtz_weights_kw_singular)(SQT_REAL *ce,
							  gint ne, gint Nk,
							  SQT_REAL *Kq,
							  gint N,
							  SQT_REAL k,
							  SQT_REAL s0,
							  SQT_REAL t0,
							  SQT_REAL *w,
							  SQT_REAL *work)

/*workspace size 3*ne + 12*ne*/

{
  /* SQT_REAL Knm[2048], x[3], n[3], J ; */
  SQT_REAL x[3], n[3], J, *swork ;
  gpointer data[SQT_DATA_WIDTH] ;

  g_assert(work != NULL) ;
  swork = &(work[3*ne]) ;
  SQT_FUNCTION_NAME(sqt_element_interp)(ce, ne, Nk, s0, t0, x, n, &J,
					NULL, work) ;
  
  data[SQT_DATA_TARGET]     = x ; 
  data[SQT_DATA_MATRIX]     = Kq ;
  data[SQT_DATA_KNM]        = work ;
  data[SQT_DATA_NKNM]       = &ne ;
  data[SQT_DATA_ORDER_K]    = &Nk ;
  data[SQT_DATA_WAVENUMBER] = &k ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)helmholtz_quad_weights ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)helmholtz_quad_weights ;
#endif /*SQT_SINGLE_PRECISION*/

  memset(w, 0, 4*ne*sizeof(SQT_REAL)) ;
  SQT_FUNCTION_NAME(sqt_singular_quad_kw_vector)(ce, ne, Nk, s0, t0, N, func,
						 w, 4*ne, data, swork) ;
  
  return 0 ;
}

#if 0
gint
SQT_FUNCTION_NAME(sqt_helmholtz_source_target_kw_self)(SQT_REAL *xe,
						       gint xstr, gint ne,
						       SQT_REAL *K,
						       gint nK,
						       gint N,
						       SQT_REAL k,
						       SQT_REAL *s,
						       gint sstr,
						       SQT_REAL *t,
						       gint tstr,
						       SQT_REAL *Ast,
						       SQT_REAL *work)

/*workspace size: 3*ne + ne*/
  
{
  gint i, i3 = 3 ;
  SQT_REAL *ce, al, bt, *swork ;

  ce = work ; swork = &(ce[3*ne]) ;
  
#ifndef SQT_SINGLE_PRECISION
  /*element geometric interpolation coefficients*/
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, ne, i3, ne, al, K, ne, xe, xstr, bt, ce, i3) ;
#else
  g_assert_not_reached() ; 
  al = 1.0 ; bt = 0.0 ;
  blaswrap_sgemm(FALSE, FALSE, ne, i3, ne, al, K, ne, xe, xstr, bt, ce, i3) ;
#endif /*SQT_SINGLE_PRECISION*/

  memset(Ast, 0, 2*ne*ne*sizeof(SQT_REAL)) ;
  for ( i = 0 ; i < ne ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_laplace_weights_kw_singular)(ce, ne, nK, K, N,
						       s[i*sstr],
						       t[i*tstr],
						       &(Ast[i*2*ne]),
						       swork) ;
  }
  
  return 0 ;
}
#endif
