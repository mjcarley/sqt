/* This file is part of SQT, a library for Singular Quadrature on Triangles
 *
 * Copyright (C) 2023 Michael Carley
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

/* #define BLOCK_INDEXED_QUADRATURE */

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
  SQT_REAL *work  = data[SQT_DATA_KNM] ;
  SQT_REAL k = *((SQT_REAL *)data[SQT_DATA_WAVENUMBER]) ;
  SQT_REAL R, dR, G[2], dG[2], E[2], wt, d1 = 1.0, d0 = 0.0 ;
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
  blaswrap_dgemv(TRUE, nq, nq, d1, Kq, nq, Knm, i1, d0, work, i1) ;
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ;
  blaswrap_sgemv(TRUE, nq, nq, d1, Kq, nq, Knm, i1, d0, work, i1) ;
#endif /*SQT_SINGLE_PRECISION*/
  
#ifndef SQT_SINGLE_PRECISION
  wt = w*G[0] ;
  blaswrap_daxpy(nq, wt, work, i1, &(quad[     0]), i2) ;
  wt = w*G[1] ;
  blaswrap_daxpy(nq, wt, work, i1, &(quad[     1]), i2) ;
  wt = w*dG[0] ;
  blaswrap_daxpy(nq, wt, work, i1, &(quad[nc/2+0]), i2) ;
  wt = w*dG[1] ;
  blaswrap_daxpy(nq, wt, work, i1, &(quad[nc/2+1]), i2) ;
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ; /*untested code*/
  /* wt = w*G ; */
  /* blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   0]), i1) ; */
  /* wt = w*dG ; */
  /* blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[nc/2]), i1) ; */
#endif /*SQT_SINGLE_PRECISION*/
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_helmholtz_weights_kw_adaptive)(SQT_REAL k,
							  SQT_REAL *ce,
							  gint ne, gint Nk,
							  SQT_REAL *Kq,
							  SQT_REAL *q, gint nq,	
							  SQT_REAL tol,
							  gint dmax,
							  SQT_REAL *x,
							  SQT_REAL *w,
							  SQT_REAL *work)

/*
 * workspace size: 16*dmax*ne + 12*ne + ne
 */
{
  gpointer data[SQT_DATA_WIDTH] ;
  
  data[SQT_DATA_TARGET]     = x ; 
  data[SQT_DATA_WAVENUMBER] = &k ; 
  data[SQT_DATA_MATRIX]     = Kq ;
  data[SQT_DATA_KNM]        = &(work[16*dmax*ne+12*ne]) ;
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

static gint helmholtz_quad_weights_tri(SQT_REAL s, SQT_REAL t,
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
  SQT_REAL k = *((SQT_REAL *)data[SQT_DATA_WAVENUMBER]) ;
  SQT_REAL R, dR, G[2], dG[2], E[2], wt, d1 = 1.0 ;
  gint i1 = 1, i2 = 2 ;

  /*Koornwinder polynomials at evaluation point*/
  SQT_FUNCTION_NAME(sqt_koornwinder_nm)(Nk, s, t, Knm, 1, nq) ;

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

gint SQT_FUNCTION_NAME(sqt_helmholtz_weights_tri_singular)(SQT_REAL k,
							   SQT_REAL *xe,
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
  data[SQT_DATA_WAVENUMBER] = &k ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)helmholtz_quad_weights_tri ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)helmholtz_quad_weights_tri ;
#endif /*SQT_SINGLE_PRECISION*/

  memset(w, 0, 4*ne*sizeof(SQT_REAL)) ;
  SQT_FUNCTION_NAME(sqt_singular_quad_tri)(xe, xstr, ne, s0, t0, N, func,
					   w, 4*nqk, data) ;
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_helmholtz_weights_tri_adaptive)(SQT_REAL k,
							   SQT_REAL *xe,
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
  data[SQT_DATA_WAVENUMBER] = &k ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)helmholtz_quad_weights_tri ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)helmholtz_quad_weights_tri ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_adaptive_quad_tri)(xe, xstr, ne, q, nq, func,
					   w, 4*nqk, tol, dmax, data) ;
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_helmholtz_weights_kw_singular)(SQT_REAL k,
							  SQT_REAL *ce,
							  gint ne, gint Nk,
							  SQT_REAL *Kq,
							  gint N,
							  SQT_REAL s0,
							  SQT_REAL t0,
							  SQT_REAL *w,
							  SQT_REAL *work)

/*workspace size 3*ne + 12*ne + ne*/

{
  SQT_REAL x[3], n[3], J, *swork ;
  gpointer data[SQT_DATA_WIDTH] ;

  g_assert(work != NULL) ;
  swork = &(work[3*ne]) ;
  SQT_FUNCTION_NAME(sqt_element_interp)(ce, ne, Nk, s0, t0, x, n, &J,
					NULL, work) ;
  
  data[SQT_DATA_TARGET]     = x ; 
  data[SQT_DATA_MATRIX]     = Kq ;
  data[SQT_DATA_KNM]        = &(swork[12*ne]) ;
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

gint
SQT_FUNCTION_NAME(sqt_helmholtz_source_target_tri_adaptive)(SQT_REAL k,
							    SQT_REAL *xse,
							    gint xsstr,
							    gint nse,
							    SQT_REAL *q,
							    gint nq,
							    SQT_REAL *Kq,
							    gint nqk, gint nK,
							    SQT_REAL tol,
							    gint dmax,
							    SQT_REAL *xte,
							    gint xtstr,
							    gint nte,
							    SQT_REAL *s,
							    gint sstr,
							    SQT_REAL *t,
							    gint tstr,
							    gint nt,
							    SQT_REAL *Ast)
{
  gint i ;
  SQT_REAL x[3], n[3], J ;

  memset(Ast, 0, 4*nqk*nt*sizeof(SQT_REAL)) ;
  for ( i = 0 ; i < nt ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_element_point_3d)(xte, xtstr, nte,
					    s[i*sstr], t[i*tstr], x, n, &J) ;
    SQT_FUNCTION_NAME(sqt_helmholtz_weights_tri_adaptive)(k, xse, xsstr, nse,
							  q, nq, Kq, nqk, nK,
							  tol, dmax, x,
							  &(Ast[i*4*nqk])) ;
  }
  
  return 0 ;
}

gint
SQT_FUNCTION_NAME(sqt_helmholtz_source_target_kw_adaptive)(SQT_REAL k,
							   SQT_REAL *xse,
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
  
  memset(Ast, 0, 4*nse*nte*sizeof(SQT_REAL)) ;  
  SQT_FUNCTION_NAME(sqt_helmholtz_matrix_kw_adaptive)(k, ce, nse, Ns, Ks,
						      q, nq,
						      tol, dmax,
						      xte, tstr, nte,
						      Ast, work) ;  
  return 0 ;
}

static gint helmholtz_quad_matrix(SQT_REAL s, SQT_REAL t,
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
  SQT_REAL *work  = data[SQT_DATA_WORK] ;
  SQT_REAL k = *((SQT_REAL *)data[SQT_DATA_WAVENUMBER]) ;
  SQT_REAL G[2], dG[2], E[2], d0 = 0.0, d1 = 1.0, R, dR ;
  SQT_REAL r[500], wt ;
  gint i1 = 1, i2 = 2, i, ns ;

#ifndef SQT_SINGLE_PRECISION  
  blaswrap_dgemv(TRUE, nq, nq, d1, Kq, nq, Knm, i1, d0, work, i1) ;
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ;
  blaswrap_sgemv(TRUE, nq, nq, d1, Kq, nq, Knm, i1, d0, work, i1) ;
#endif /*SQT_SINGLE_PRECISION*/
  
  ns = nc/nx ;
  for ( i = 0 ; i < nx ; i ++ ) {
    sqt_vector_diff(&(r[3*i]), &(x[i*xstr]), y) ;
    R = sqrt(r[3*i+0]*r[3*i+0] + r[3*i+1]*r[3*i+1] + r[3*i+2]*r[3*i+2]) ;
    dR = sqt_vector_diff_scalar(&(x[i*xstr]), y, n)/R ;
    E[0] = cos(k*R) ; E[1] = sin(k*R) ;
    G[0] = 0.25*M_1_PI/R ;
    G[1] = G[0]*E[1] ;

    dG[0] = (E[0] + k*R*E[1])*G[0]*dR/R ;
    dG[1] = (E[1] - k*R*E[0])*G[0]*dR/R ;

    G[0] *= E[0] ;
#ifndef SQT_SINGLE_PRECISION
  /*it should be possible to do this with one matrix multiplication*/
  wt = w*G[0] ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1,
		 &(quad[i*ns +  0]), i2) ;
  wt = w*G[1] ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1,
		 &(quad[i*ns +  1]), i2) ;
  /*and this*/
  wt = w*dG[0] ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1,
		 &(quad[i*ns + ns/2 +  0]), i2) ;
  wt = w*dG[1] ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1,
		 &(quad[i*ns + ns/2 +  1]), i2) ;

#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ; /*untested code*/
  /* wt = w*G ; */
  /* blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   0]), i1) ; */
  /* wt = w*dG ; */
  /* blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[nc/2]), i1) ; */
#endif
  }

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_helmholtz_matrix_kw_adaptive)(SQT_REAL k,
							 SQT_REAL *ce,
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
  data[SQT_DATA_WORK]      = &(work[4*dmax*2*ne*nx]) ;
  data[SQT_DATA_WAVENUMBER] = &k ;
				    
/* #define USE_VECTOR */
#ifdef SQT_SINGLE_PRECISION
#ifdef USE_VECTOR
  sqt_quadrature_vec_func_f_t func =
    (sqt_quadrature_vec_func_f_t)laplace_quad_vec_matrix ;
#else /*USE_VECTOR*/
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)helmholtz_quad_matrix ;
#endif /*USE_VECTOR*/
#else /*SQT_SINGLE_PRECISION*/
#ifdef USE_VECTOR
  sqt_quadrature_vec_func_t func =
    (sqt_quadrature_vec_func_t)laplace_quad_vec_matrix ;
#else /*USE_VECTOR*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)helmholtz_quad_matrix ;
#endif /*USE_VECTOR*/
#endif /*SQT_SINGLE_PRECISION*/

#ifdef USE_VECTOR
  SQT_FUNCTION_NAME(sqt_adaptive_quad_vec_kw)(ce, ne, Nk, q, nq, func,
  					      Ast, 4*ne*nx, tol, dmax,
  					      data, work) ;
#else /*USE_VECTOR*/
  /* SQT_FUNCTION_NAME(sqt_adaptive_quad_kw)(ce, ne, Nk, q, nq, func, */
  /* 					  Ast, 2*ne*nx, tol, dmax, */
  /* 					  data, work) ; */
  SQT_FUNCTION_NAME(sqt_tree_quad_kw)(ce, ne, Nk, q, nq, func,
				      Ast, 4*ne*nx, tol, dmax,
				      data, work) ;
#endif /*USE_VECTOR*/
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_helmholtz_weights_tri_basic)(SQT_REAL k,
							SQT_REAL *xe,
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
  data[SQT_DATA_WAVENUMBER] = &k ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)helmholtz_quad_weights_tri ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)helmholtz_quad_weights_tri ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_basic_quad_tri)(xe, xstr, ne, q, nq, func,
					w, 4*nqk, data) ;
  
  return 0 ;
}

gint
SQT_FUNCTION_NAME(sqt_helmholtz_source_target_tri_basic)(SQT_REAL k,
							 SQT_REAL *xse,
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

  memset(Ast, 0, 4*nqk*nt*sizeof(SQT_REAL)) ;
  for ( i = 0 ; i < nt ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_element_point_3d)(xte, xtstr, nte,
					    s[i*sstr], t[i*tstr], x, n, &J) ;
    SQT_FUNCTION_NAME(sqt_helmholtz_weights_tri_basic)(k, xse, xsstr, nse,
						       q, nq, Kq, nqk, nK,
						       x, &(Ast[i*4*nqk])) ;
  }
  
  return 0 ;
}

#ifndef BLOCK_INDEXED_QUADRATURE
static gint helmholtz_quad_matrix_indexed(SQT_REAL s, SQT_REAL t,
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
  SQT_REAL *work  = data[SQT_DATA_KNM] ;
  SQT_REAL k = *((SQT_REAL *)data[SQT_DATA_WAVENUMBER]) ;
  gint *idx =           data[SQT_DATA_INDICES] ;
  SQT_REAL G[2], dG[2], E[2], d0 = 0.0, d1 = 1.0, R, dR ;
  SQT_REAL r[3], wt ;
  gint i1 = 1, i2 = 2, i, j, ns ;

#ifndef SQT_SINGLE_PRECISION
  blaswrap_dgemv(TRUE, nq, nq, d1, Kq, nq, Knm, i1, d0, work, i1) ;
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ;
  blaswrap_sgemv(TRUE, nq, nq, d1, Kq, nq, Knm, i1, d0, work, i1) ;
#endif /*SQT_SINGLE_PRECISION*/
  
  ns = nc/nx ;
  for ( j = 0 ; j < nx ; j ++ ) {
    i = idx[j] ;
    sqt_vector_diff(r, &(x[i*xstr]), y) ;
    /* R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) ; */
    /* dR = sqt_vector_diff_scalar(&(x[i*xstr]), y, n)/R ; */
    R = r[0]*r[0] + r[1]*r[1] + r[2]*r[2] ;
    /* dR = sqt_vector_diff_scalar(&(x[i*xstr]), y, n)/R ; */
    dR = sqt_vector_scalar(r, n)/R ;
    R = sqrt(R) ;
    E[0] = cos(k*R) ; E[1] = sin(k*R) ;
    G[0] = 0.25*M_1_PI/R ;
    G[1] = G[0]*E[1] ;

    /* dG[0] = (E[0] + k*R*E[1])*G[0]*dR/R ; */
    /* dG[1] = (E[1] - k*R*E[0])*G[0]*dR/R ; */
    dG[0] = (E[0] + k*R*E[1])*G[0]*dR ;
    dG[1] = (E[1] - k*R*E[0])*G[0]*dR ;

    G[0] *= E[0] ;
#ifndef SQT_SINGLE_PRECISION
    wt = w*G[0] ;
    blaswrap_daxpy(nq, wt, work, i1, &(quad[j*ns+0]), i2) ;
    wt = w*G[1] ;
    blaswrap_daxpy(nq, wt, work, i1, &(quad[j*ns+1]), i2) ;
    wt = w*dG[0] ;
    blaswrap_daxpy(nq, wt, work, i1, &(quad[j*ns+ns/2+0]), i2) ;
    wt = w*dG[1] ;
    blaswrap_daxpy(nq, wt, work, i1, &(quad[j*ns+ns/2+1]), i2) ;
  /*
   * this is significantly (factor of four) slower than four vector
   * multiplications
   */
  /* blaswrap_dgemm(FALSE, FALSE, nq, i2, i1, w, work, i1, G, i2, d1, */
  /* 		 &(quad[j*ns]), i2) ; */
  /* blaswrap_dgemm(FALSE, FALSE, nq, i2, i1, w, work, i1, dG, i2, d1, */
  /* 		 &(quad[j*ns+ns/2]), i2) ; */
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ; /*untested code*/
  /* wt = w*G ; */
  /* blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   0]), i1) ; */
  /* wt = w*dG ; */
  /* blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[nc/2]), i1) ; */
#endif
  }
  
  return 0 ;
}

#else /*BLOCK_INDEXED_QUADRATURE*/

static gint helmholtz_quad_matrix_block(SQT_REAL s, SQT_REAL t,
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
  SQT_REAL *work  = data[SQT_DATA_KNM] ;
  SQT_REAL k = *((SQT_REAL *)data[SQT_DATA_WAVENUMBER]) ;
  SQT_REAL G[2], dG[2], E[2], d0 = 0.0, d1 = 1.0, R, dR ;
  SQT_REAL r[3], wt ;
  gint i1 = 1, i2 = 2, i, j, ns ;

#ifndef SQT_SINGLE_PRECISION
  blaswrap_dgemv(TRUE, nq, nq, d1, Kq, nq, Knm, i1, d0, work, i1) ;
#else /*SQT_SINGLE_PRECISION*/
  g_assert_not_reached() ;
  blaswrap_sgemv(TRUE, nq, nq, d1, Kq, nq, Knm, i1, d0, work, i1) ;
#endif /*SQT_SINGLE_PRECISION*/
  
  ns = nc/nx ;
  for ( i = 0 ; i < nx ; i ++ ) {
    sqt_vector_diff(r, &(x[i*xstr]), y) ;
    R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) ;
    dR = sqt_vector_diff_scalar(&(x[i*xstr]), y, n)/R ;
    E[0] = cos(k*R) ; E[1] = sin(k*R) ;
    G[0] = 0.25*M_1_PI/R ;
    G[1] = G[0]*E[1] ;

    dG[0] = (E[0] + k*R*E[1])*G[0]*dR/R ;
    dG[1] = (E[1] - k*R*E[0])*G[0]*dR/R ;

    G[0] *= E[0] ;
#ifndef SQT_SINGLE_PRECISION
    wt = w*G[0] ;
    blaswrap_daxpy(nq, wt, work, i1, &(quad[i*ns+0]), i2) ;
    wt = w*G[1] ;
    blaswrap_daxpy(nq, wt, work, i1, &(quad[i*ns+1]), i2) ;
    wt = w*dG[0] ;
    blaswrap_daxpy(nq, wt, work, i1, &(quad[i*ns+ns/2+0]), i2) ;
    wt = w*dG[1] ;
    blaswrap_daxpy(nq, wt, work, i1, &(quad[i*ns+ns/2+1]), i2) ;
#else /*SQT_SINGLE_PRECISION*/
    g_assert_not_reached() ; /*untested code*/
    /* wt = w*G ; */
    /* blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[   0]), i1) ; */
    /* wt = w*dG ; */
    /* blaswrap_sgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, &(quad[nc/2]), i1) ; */
#endif
  }
  
  return 0 ;
}
#endif /*BLOCK_INDEXED_QUADRATURE*/

gint
SQT_FUNCTION_NAME(sqt_helmholtz_source_indexed_kw_adaptive)(SQT_REAL k,
							    SQT_REAL *xse,
							    gint sstr, gint nse,
							    SQT_REAL *q,
							    gint nq,
							    SQT_REAL *Ks,
							    gint Ns,
							    SQT_REAL tol,
							    gint dmax,
							    SQT_REAL *xte,
							    gint tstr,
							    gint *idx,
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

  data[SQT_DATA_NUMBER]    = &nte ;
  data[SQT_DATA_MATRIX]    = Ks ;
  data[SQT_DATA_KNM]       = Knm ;
  data[SQT_DATA_NKNM]      = &nse ;
  data[SQT_DATA_ORDER_K]   = &Ns ;
  data[SQT_DATA_N_QUAD_POINTS] = &nq ;
  data[SQT_DATA_WAVENUMBER] = &k ;

#ifdef BLOCK_INDEXED_QUADRATURE
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)helmholtz_quad_matrix_block ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)helmholtz_quad_matrix_block ;
#endif /*SQT_SINGLE_PRECISION*/

/*
 * work space size: 3*nse + (12*nse + 2*nte) + (4*dmax*4*nse*nte+12*nse) +
 * str*nte + 4*nse*nte
 */
  
  str = 3 ;
  memset(Ast, 0, 4*nse*nte*sizeof(SQT_REAL)) ;
  xb = &(awork[4*dmax*4*nse*nte+12*nse]) ;
  A = &(xb[str*nte]) ;
  memset(A, 0, 4*nse*nte*sizeof(SQT_REAL)) ;

  for ( i = 0 ; i < nte ; i ++ ) {
    j = idx[i] ; 
    xb[str*i+0] = xte[tstr*j+0] ; 
    xb[str*i+1] = xte[tstr*j+1] ; 
    xb[str*i+2] = xte[tstr*j+2] ; 
  }

  data[SQT_DATA_TARGET]    = xb ; 
  data[SQT_DATA_STRIDE]    = &str ;   
  
  SQT_FUNCTION_NAME(sqt_tree_quad_kw)(ce, nse, Ns, q, nq, func,
  				      A, 4*nse*nte, tol, dmax,
  				      data, awork) ;

  for ( i = 0 ; i < nte ; i ++ ) {
    j = idx[i] ; 
    for ( n = 0 ; n < 4*nse ; n ++ ) {
      Ast[j*4*nse + n] = A[i*4*nse + n] ;
    }
  }

#else /*BLOCK_INDEXED_QUADRATURE*/
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)helmholtz_quad_matrix_indexed ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)helmholtz_quad_matrix_indexed ;
#endif /*SQT_SINGLE_PRECISION*/

  data[SQT_DATA_TARGET]    = xte ; 
  data[SQT_DATA_STRIDE]    = &tstr ; 
  data[SQT_DATA_INDICES]   = idx ;

  memset(Ast, 0, 4*nse*nte*sizeof(SQT_REAL)) ;
  SQT_FUNCTION_NAME(sqt_tree_quad_kw)(ce, nse, Ns, q, nq, func,
  				      Ast, 4*nse*nte, tol, dmax,
  				      data, awork) ;
#endif /*BLOCK_INDEXED_QUADRATURE*/
  
return 0 ;
}

gint
SQT_FUNCTION_NAME(sqt_helmholtz_source_target_kw_self)(SQT_REAL k,
						       SQT_REAL *xe,
						       gint xstr, gint ne,
						       SQT_REAL *K,
						       gint nK,
						       gint N,
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

  memset(Ast, 0, 4*ne*ne*sizeof(SQT_REAL)) ;
  for ( i = 0 ; i < ne ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_helmholtz_weights_kw_singular)(k,
							 ce, ne, nK, K, N,
							 s[i*sstr],
							 t[i*tstr],
							 &(Ast[i*4*ne]),
							 swork) ;
  }
  
  return 0 ;
}


