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
