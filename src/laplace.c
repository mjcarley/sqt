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

static gint SQT_FUNCTION_NAME(laplace_quad_weights)(SQT_REAL s, SQT_REAL t,
						    SQT_REAL w,
						    SQT_REAL *y, SQT_REAL *n,
						    SQT_REAL *quad, gint nc,
						    gpointer data[])

{
  SQT_REAL *Knm = data[SQT_DATA_KNM] ;
  gint nq = *((gint *)(data[SQT_DATA_NKNM])) ;
  gint Nk = *((gint *)(data[SQT_DATA_ORDER_K])) ;
  SQT_REAL *x  = data[SQT_DATA_TARGET] ;
  SQT_REAL R, G, dG ;
  SQT_REAL *Kq  = data[SQT_DATA_MATRIX] ;
  SQT_REAL wt, d1 = 1.0 ;
  gint i1 = 1 ;

  /*Koornwinder polynomials at evaluation point*/
  SQT_FUNCTION_NAME(sqt_koornwinder_nm)(Nk, s, t, 1, nq, Knm) ;
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
#endif
  
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
    (sqt_quadrature_func_f_t)laplace_quad_weights_f ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_adaptive_quad_tri)(xe, xstr, ne, q, nq, func,
					   w, 2*nqk, tol, dmax, data) ;
  
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
  SQT_REAL Knm[2048], x[3], n[3], J ;
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
    (sqt_quadrature_func_f_t)laplace_quad_weights_f ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_singular_quad_tri)(xe, xstr, ne, s0, t0, N, func,
					   w, 2*nqk, data) ;
  
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
SQT_FUNCTION_NAME(sqt_laplace_source_target_tri_self)(SQT_REAL *xe,
						      gint xstr, gint ne,
						      SQT_REAL *Kq,
						      gint nqk, gint nK,
						      gint N,
						      SQT_REAL *s,
						      gint sstr,
						      SQT_REAL *t,
						      gint tstr,
						      SQT_REAL *Ast)
{
  gint i ;

  memset(Ast, 0, 2*nqk*nqk*sizeof(SQT_REAL)) ;
  for ( i = 0 ; i < nqk ; i ++ ) {
    /* fprintf(stderr, "%d %lg %lg %lg\n", */
    /* 	    i, s[i*sstr], t[i*tstr], s[i*sstr] + t[i*tstr] - 1.0) ; */
    SQT_FUNCTION_NAME(sqt_laplace_weights_tri_singular)(xe, xstr, ne,
							Kq, nqk, nK, N,
							s[i*sstr],
							t[i*tstr],
							&(Ast[i*2*nqk])) ;
  }
  
  return 0 ;
}
