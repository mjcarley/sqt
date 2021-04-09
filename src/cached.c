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

#include "config.h"

static void cached_quad_kw(SQT_REAL *ce, gint ne, gint Nk,
			   SQT_REAL *st, SQT_REAL wt,
			   SQT_REAL *q, gint nq, 
#ifdef SQT_SINGLE_PRECISION
			   sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
			   sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
			   SQT_REAL *quad, gint nc,
			   gint *icache, SQT_REAL *xcache, gint cstr,
			   gint d, gint b,
			   gpointer data, SQT_REAL *work)

{
  gint i, off, init ;
  /* SQT_REAL s, t, J, *y, *n, *xc, si, ti, work[3*453], w ; */
  SQT_REAL s, t, J, *y, *n, *xc, si, ti, w ;

  /* cstr = NBI_CACHE_STRIDE + nq ; */
  off = sqt_cache_level_offset(d) ;
  xc = &(xcache[(off+b)*cstr*nq]) ;
  init = 0 ;
  
  if ( icache[off+b] == 0 ) {
    init = 1 ;
    for ( i = 0 ; i < nq ; i ++ ) {
      si = q[3*i+0] ; ti = q[3*i+1] ;
      /*interpolating coordinates on unit triangle*/
      s = sqt_cache_coordinate_s(xc,cstr,i) =
	(1.0-si-ti)*st[0] + si*st[2] + ti*st[4] ;
      t = sqt_cache_coordinate_t(xc,cstr,i) =
	(1.0-si-ti)*st[1] + si*st[3] + ti*st[5] ;
      y = sqt_cache_position(xc,cstr,i) ;
      n = sqt_cache_normal(xc,cstr,i) ;

      /*coordinates on physical triangle*/
      SQT_FUNCTION_NAME(sqt_element_interp)(ce, ne, Nk, s, t, y, n, &J,
					    NULL, work) ;
      sqt_cache_weight(xc,cstr,i) = q[3*i+2]*J*wt ;
    }
    icache[off+b] = 1 ;
  }

  /* g_assert(icache[off+b] == 1) ; */
  
  for ( i = 0 ; i < nq ; i ++ ) {
    s = sqt_cache_coordinate_s(xc,cstr,i) ;
    t = sqt_cache_coordinate_t(xc,cstr,i) ;
    w = sqt_cache_weight(xc,cstr,i) ;
    y = sqt_cache_position(xc,cstr,i) ;
    n = sqt_cache_normal(xc,cstr,i) ;

    func(s, t, w, y, n, NULL, 0, quad, nc, init, data) ;
    /* if ( isnan(quad[0]) ) */
    /*   g_error("NaN") ; */
  }

  return ;
}

static gint cached_quad_kw_recursion(SQT_REAL *ce, gint ne, gint Nk,
				     SQT_REAL *st, SQT_REAL wt,
				     SQT_REAL *q, gint nq,
#ifdef SQT_SINGLE_PRECISION
				     sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
				     sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
				     SQT_REAL *quad, gint nc,
				     SQT_REAL tol,
				     gint d, gint dmax, gint b,
				     gint *icache, SQT_REAL *xcache, gint cstr,
				     gpointer data, SQT_REAL *work)

{
  gint i ;
  SQT_REAL *q0, *q1, *q2, *q3 ;
  SQT_REAL st0[6], st1[6], st2[6], st3[6], *wtmp ;
  gboolean recurse ;

  if ( d == dmax ) return 0 ;

  memset(work, 0, 4*nc*sizeof(SQT_REAL)) ;
  q0 = &(work[0]) ; q1 = &(q0[nc]) ; q2 = &(q1[nc]) ; q3 = &(q2[nc]) ;
  wtmp = &(q3[nc]) ;
  
  wt *= 0.25 ;
  
  sqt_triangle_divide_loop30(st, st0) ;
  cached_quad_kw(ce, ne, Nk, st0, wt, q, nq, func, q0, nc, icache, xcache,
		 cstr, d+1, 4*b+0, data, wtmp) ;
  sqt_triangle_divide_loop31(st, st1) ;
  cached_quad_kw(ce, ne, Nk, st1, wt, q, nq, func, q1, nc, icache, xcache,
		 cstr, d+1, 4*b+1, data, wtmp) ;
  sqt_triangle_divide_loop32(st, st2) ;
  cached_quad_kw(ce, ne, Nk, st2, wt, q, nq, func, q2, nc, icache, xcache,
		 cstr, d+1, 4*b+2, data, wtmp) ;
  sqt_triangle_divide_loop33(st, st3) ;
  cached_quad_kw(ce, ne, Nk, st3, wt, q, nq, func, q3, nc, icache, xcache,
		 cstr, d+1, 4*b+3, data, wtmp) ;

  recurse = FALSE ;

  for ( i = 0 ; i < nc ; i ++ ) {
    if ( fabs(quad[i] - q0[i] - q1[i] - q2[i] - q3[i]) > tol ) {
      recurse = TRUE ; break ;
    }
  }

  if ( !recurse ) return 0 ;

  cached_quad_kw_recursion(ce, ne, Nk, st0, wt, q, nq, func,
			   q0, nc, tol, d+1, dmax, 4*b+0, icache, xcache,
			   cstr, data, &(work[4*nc])) ;
  cached_quad_kw_recursion(ce, ne, Nk, st1, wt, q, nq, func,
			   q1, nc, tol, d+1, dmax, 4*b+1, icache, xcache,
			   cstr, data, &(work[4*nc])) ;
  cached_quad_kw_recursion(ce, ne, Nk, st2, wt, q, nq, func,
			   q2, nc, tol, d+1, dmax, 4*b+2, icache, xcache,
			   cstr, data, &(work[4*nc])) ;
  cached_quad_kw_recursion(ce, ne, Nk, st3, wt, q, nq, func,
			   q3, nc, tol, d+1, dmax, 4*b+3, icache, xcache,
			   cstr, data, &(work[4*nc])) ;
  for ( i = 0 ; i < nc ; i ++ ) {
    quad[i] = q0[i] + q1[i] + q2[i] + q3[i] ;
  }

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_cached_quad_kw)(SQT_REAL *ce, gint ne, gint Nk,
					   SQT_REAL *q, gint nq,
#ifdef SQT_SINGLE_PRECISION
					   sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
					   sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
					   SQT_REAL *quad, gint nc,
					   SQT_REAL tol, gint dmax,
					   gint *icache, SQT_REAL *xcache,
					   gint cstr,
					   gboolean init,
					   gpointer data, SQT_REAL *work)

/*
 * work space size: 4*dmax*nc + 3*nq
 */

{
  SQT_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0} ;

  g_assert(cstr >= SQT_CACHE_STRIDE) ;
  
  if ( init ) memset(icache, 0, sqt_cache_level_offset(dmax+1)*sizeof(gint)) ;
  
  memset(quad, 0, nc*sizeof(SQT_REAL)) ;
  cached_quad_kw(ce, ne, Nk, st, 1.0, q, nq, func, quad, nc,
		 icache, xcache, cstr, 0, 0, data, work) ;

  return cached_quad_kw_recursion(ce, ne, Nk, st, 1.0, q, nq,
				  func, quad, nc, tol, 0, dmax, 0,
				  icache, xcache, cstr,
				  data, work) ;
}
