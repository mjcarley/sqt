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

static gint laplace_quad_weights(SQT_REAL s, SQT_REAL t,
				 SQT_REAL w,
				 SQT_REAL *y, SQT_REAL *n,
				 SQT_REAL *quad, gint nc,
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
  /*this is to keep the compiler quiet about unused variables*/
  wt = 0.0 ; Kq[0] = sin(wt*d1*i1*dG) ;
  g_assert_not_reached() ;
#endif
  
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
    (sqt_quadrature_func_f_t)laplace_quad_weights ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights ;
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
    (sqt_quadrature_func_f_t)laplace_quad_weights ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_adaptive_quad_tri)(xe, xstr, ne, q, nq, func,
					   w, 2*nqk, tol, dmax, data) ;
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_laplace_weights_kw_adaptive)(SQT_REAL *ce,
							gint ne, gint Nk,
							SQT_REAL *Kq,
							SQT_REAL *q, gint nq,
							SQT_REAL tol,
							gint dmax,
							SQT_REAL *x,
							SQT_REAL *w)

{
  SQT_REAL Knm[2048] ;
  gpointer data[SQT_DATA_WIDTH] ;
  
  data[SQT_DATA_TARGET]    = x ; 
  data[SQT_DATA_MATRIX]    = Kq ;
  data[SQT_DATA_KNM]       = Knm ;
  data[SQT_DATA_NKNM]      = &ne ;
  data[SQT_DATA_ORDER_K]   = &Nk ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)laplace_quad_weights ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_adaptive_quad_kw)(ce, ne, Nk, q, nq, func,
					  w, 2*ne, tol, dmax, data) ;
  
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
    (sqt_quadrature_func_f_t)laplace_quad_weights ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_singular_quad_tri)(xe, xstr, ne, s0, t0, N, func,
					   w, 2*nqk, data) ;
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_laplace_weights_kw_singular)(SQT_REAL *ce,
							gint ne, gint Nk,
							SQT_REAL *Kq,
							gint N,
							SQT_REAL s0,
							SQT_REAL t0,
							SQT_REAL *w)

{
  SQT_REAL Knm[2048], x[3], n[3], J ;
  gpointer data[SQT_DATA_WIDTH] ;

  SQT_FUNCTION_NAME(sqt_element_interp)(ce, ne, Nk, s0, t0, x, n, &J,
					NULL, Knm) ;
  
  data[SQT_DATA_TARGET]    = x ; 
  data[SQT_DATA_MATRIX]    = Kq ;
  data[SQT_DATA_KNM]       = Knm ;
  data[SQT_DATA_NKNM]      = &ne ;
  data[SQT_DATA_ORDER_K]   = &Nk ;
  
#ifdef SQT_SINGLE_PRECISION
  sqt_quadrature_func_f_t func =
    (sqt_quadrature_func_f_t)laplace_quad_weights ;
#else /*SQT_SINGLE_PRECISION*/
  sqt_quadrature_func_t func =
    (sqt_quadrature_func_t)laplace_quad_weights ;
#endif /*SQT_SINGLE_PRECISION*/

  SQT_FUNCTION_NAME(sqt_singular_quad_kw)(ce, ne, Nk, s0, t0, N, func,
					  w, 2*ne, data) ;
  
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
							 SQT_REAL *Ast)

{
  gint i, i3 = 3 ;
  SQT_REAL ce[3*453], al, bt ;

#ifndef SQT_SINGLE_PRECISION
  /*element geometric interpolation coefficients*/
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, nse, i3, nse, al, Ks, nse,
		 xse, sstr, bt, ce, i3) ;
#else
  g_assert_not_reached() ; 
#endif /*SQT_SINGLE_PRECISION*/
  
  memset(Ast, 0, 2*nse*nte*sizeof(SQT_REAL)) ;
  for ( i = 0 ; i < nte ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_laplace_weights_kw_adaptive)(ce, nse, Ns, Ks,
						       q, nq,
						       tol, dmax,
						       &(xte[i*tstr]),
						       &(Ast[i*2*nse])) ;
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

gint
SQT_FUNCTION_NAME(sqt_laplace_source_target_kw_self)(SQT_REAL *xe,
						     gint xstr, gint ne,
						     SQT_REAL *K,
						     gint nK,
						     gint N,
						     SQT_REAL *s,
						     gint sstr,
						     SQT_REAL *t,
						     gint tstr,
						     SQT_REAL *Ast)

{
  gint i, i3 = 3 ;
  SQT_REAL ce[3*453], al, bt ;

#ifndef SQT_SINGLE_PRECISION
  /*element geometric interpolation coefficients*/
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, ne, i3, ne, al, K, ne, xe, xstr, bt, ce, i3) ;
#else
  g_assert_not_reached() ; 
#endif /*SQT_SINGLE_PRECISION*/

  memset(Ast, 0, 2*ne*ne*sizeof(SQT_REAL)) ;
  for ( i = 0 ; i < ne ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_laplace_weights_kw_singular)(ce, ne, nK, K, N,
						       s[i*sstr],
						       t[i*tstr],
						       &(Ast[i*2*ne])) ;
  }
  
  return 0 ;
}
