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

/* #define TRIANGLE_TRACE */

static void adaptive_quad_tri(SQT_REAL *xe, gint xstr, gint ne,
			      SQT_REAL *st, SQT_REAL wt,
			      SQT_REAL *q, gint nq, 
#ifdef SQT_SINGLE_PRECISION
			      sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
			      sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
			      SQT_REAL *quad, gint nc,
			      gpointer data)

{
  gint i ;
  SQT_REAL s, t, w, J, y[3], n[3], si, ti ;
  
  for ( i = 0 ; i < nq ; i ++ ) {
    si = q[3*i+0] ; ti = q[3*i+1] ; w = q[3*i+2] ;
    /*interpolating coordinates on unit triangle*/
    s = (1.0-si-ti)*st[0] + si*st[2] + ti*st[4] ;
    t = (1.0-si-ti)*st[1] + si*st[3] + ti*st[5] ;
    /*coordinates on physical triangle*/
    SQT_FUNCTION_NAME(sqt_element_point_3d)(xe, xstr, ne, s, t, y, n, &J) ;
    w *= J*wt ;
    func(s, t, w, y, n, quad, nc, data) ;
  }

  return ;
}

static gint adaptive_quad_tri_recursion(SQT_REAL *xe, gint xstr, gint ne,
					SQT_REAL *st, SQT_REAL wt,
					SQT_REAL *q, gint nq,
#ifdef SQT_SINGLE_PRECISION
					sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
					sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
					SQT_REAL *quad, gint nc,
					SQT_REAL tol, gint dmax,
					gpointer data)

{
  gint i ;
  SQT_REAL work[2048], *q0, *q1, *q2, *q3 ;
  SQT_REAL st0[6], st1[6], st2[6], st3[6] ;
  gboolean recurse ;

  if ( dmax == 0 ) return 0 ;

  g_assert(4*nc <= 2048) ;
  
  memset(work, 0, 4*nc*sizeof(SQT_REAL)) ;
  q0 = &(work[0]) ; q1 = &(q0[nc]) ; q2 = &(q1[nc]) ; q3 = &(q2[nc]) ;

  wt *= 0.25 ;
  
  sqt_triangle_divide_loop30(st, st0) ;
  adaptive_quad_tri(xe, xstr, ne, st0, wt, q, nq, func, q0, nc, data) ;
  sqt_triangle_divide_loop31(st, st1) ;
  adaptive_quad_tri(xe, xstr, ne, st1, wt, q, nq, func, q1, nc, data) ;
  sqt_triangle_divide_loop32(st, st2) ;
  adaptive_quad_tri(xe, xstr, ne, st2, wt, q, nq, func, q2, nc, data) ;
  sqt_triangle_divide_loop33(st, st3) ;
  adaptive_quad_tri(xe, xstr, ne, st3, wt, q, nq, func, q3, nc, data) ;

  recurse = FALSE ;

  for ( i = 0 ; i < nc ; i ++ ) {
    if ( fabs(quad[i] - q0[i] - q1[i] - q2[i] - q3[i]) > tol ) {
      recurse = TRUE ; break ;
    }
  }

  if ( !recurse ) return 0 ;

  adaptive_quad_tri_recursion(xe, xstr, ne, st0, wt, q, nq, func,
			      q0, nc, tol, dmax-1, data) ;
  adaptive_quad_tri_recursion(xe, xstr, ne, st1, wt, q, nq, func,
			      q1, nc, tol, dmax-1, data) ;
  adaptive_quad_tri_recursion(xe, xstr, ne, st2, wt, q, nq, func,
			      q2, nc, tol, dmax-1, data) ;
  adaptive_quad_tri_recursion(xe, xstr, ne, st3, wt, q, nq, func,
			      q3, nc, tol, dmax-1, data) ;
  
  for ( i = 0 ; i < nc ; i ++ ) {
    quad[i] = q0[i] + q1[i] + q2[i] + q3[i] ;
  }

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_adaptive_quad_tri)(SQT_REAL *xe, gint xstr, gint ne,
					      SQT_REAL *q, gint nq,
#ifdef SQT_SINGLE_PRECISION
					      sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
					      sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
					      SQT_REAL *quad, gint nc,
					      SQT_REAL tol, gint dmax,
					      gpointer data)

{
  SQT_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0} ;

  memset(quad, 0, nc*sizeof(SQT_REAL)) ;
  adaptive_quad_tri(xe, xstr, ne, st, 1.0, q, nq, func, quad, nc, data) ;

  return adaptive_quad_tri_recursion(xe, xstr, ne, st, 1.0, q, nq,
				     func, quad, nc, tol, dmax,
				     data) ;
}


static void adaptive_quad_kw(SQT_REAL *ce, gint ne, gint Nk,
			      SQT_REAL *st, SQT_REAL wt,
			      SQT_REAL *q, gint nq, 
#ifdef SQT_SINGLE_PRECISION
			      sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
			      sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
			      SQT_REAL *quad, gint nc,
			      gpointer data)

{
  gint i ;
  SQT_REAL s, t, w, J, y[3], n[3], si, ti, work[3*(453+1)*(453+2)/2] ;
  
  for ( i = 0 ; i < nq ; i ++ ) {
    si = q[3*i+0] ; ti = q[3*i+1] ; w = q[3*i+2] ;
    /*interpolating coordinates on unit triangle*/
    s = (1.0-si-ti)*st[0] + si*st[2] + ti*st[4] ;
    t = (1.0-si-ti)*st[1] + si*st[3] + ti*st[5] ;
    /*coordinates on physical triangle*/
    SQT_FUNCTION_NAME(sqt_element_interp)(ce, ne, Nk, s, t, y, n, &J,
					  NULL, work) ;
    w *= J*wt ;
    func(s, t, w, y, n, quad, nc, data) ;
  }

  return ;
}

static gint adaptive_quad_kw_recursion(SQT_REAL *ce, gint ne, gint Nk,
				       SQT_REAL *st, SQT_REAL wt,
				       SQT_REAL *q, gint nq,
#ifdef SQT_SINGLE_PRECISION
				       sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
				       sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
				       SQT_REAL *quad, gint nc,
				       SQT_REAL tol, gint dmax,
				       gpointer data)

{
  gint i ;
  SQT_REAL work[2048], *q0, *q1, *q2, *q3 ;
  SQT_REAL st0[6], st1[6], st2[6], st3[6] ;
  gboolean recurse ;

  if ( dmax == 0 ) return 0 ;

  g_assert(4*nc <= 2048) ;
  
  memset(work, 0, 4*nc*sizeof(SQT_REAL)) ;
  q0 = &(work[0]) ; q1 = &(q0[nc]) ; q2 = &(q1[nc]) ; q3 = &(q2[nc]) ;

  wt *= 0.25 ;
  
  sqt_triangle_divide_loop30(st, st0) ;
  adaptive_quad_kw(ce, ne, Nk, st0, wt, q, nq, func, q0, nc, data) ;
  sqt_triangle_divide_loop31(st, st1) ;
  adaptive_quad_kw(ce, ne, Nk, st1, wt, q, nq, func, q1, nc, data) ;
  sqt_triangle_divide_loop32(st, st2) ;
  adaptive_quad_kw(ce, ne, Nk, st2, wt, q, nq, func, q2, nc, data) ;
  sqt_triangle_divide_loop33(st, st3) ;
  adaptive_quad_kw(ce, ne, Nk, st3, wt, q, nq, func, q3, nc, data) ;

  recurse = FALSE ;

  for ( i = 0 ; i < nc ; i ++ ) {
    if ( fabs(quad[i] - q0[i] - q1[i] - q2[i] - q3[i]) > tol ) {
      recurse = TRUE ; break ;
    }
  }

  if ( !recurse ) return 0 ;

  adaptive_quad_kw_recursion(ce, ne, Nk, st0, wt, q, nq, func,
			     q0, nc, tol, dmax-1, data) ;
  adaptive_quad_kw_recursion(ce, ne, Nk, st1, wt, q, nq, func,
			     q1, nc, tol, dmax-1, data) ;
  adaptive_quad_kw_recursion(ce, ne, Nk, st2, wt, q, nq, func,
			     q2, nc, tol, dmax-1, data) ;
  adaptive_quad_kw_recursion(ce, ne, Nk, st3, wt, q, nq, func,
			     q3, nc, tol, dmax-1, data) ;
  
  for ( i = 0 ; i < nc ; i ++ ) {
    quad[i] = q0[i] + q1[i] + q2[i] + q3[i] ;
  }

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_adaptive_quad_kw)(SQT_REAL *ce, gint ne, gint Nk,
					     SQT_REAL *q, gint nq,
#ifdef SQT_SINGLE_PRECISION
					     sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
					     sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
					     SQT_REAL *quad, gint nc,
					     SQT_REAL tol, gint dmax,
					     gpointer data)

{
  SQT_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0} ;

  memset(quad, 0, nc*sizeof(SQT_REAL)) ;
  adaptive_quad_kw(ce, ne, Nk, st, 1.0, q, nq, func, quad, nc, data) ;

  return adaptive_quad_kw_recursion(ce, ne, Nk, st, 1.0, q, nq,
  				     func, quad, nc, tol, dmax,
  				     data) ;
}
