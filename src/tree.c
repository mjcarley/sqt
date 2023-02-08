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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <glib.h>

#include <blaswrap.h>

#include <sqt.h>

#include "sqt-private.h"

#include "config.h"

#define SQT_TRAVERSE_STACK_SIZE 64

#define sqt_tree_level_offset(_d) ((sqt_fourpow((_d))-1)/3)
#define sqt_intlog2(_n) (31 - __builtin_clz((_n)))
#define sqt_tree_level_depth(_i) (sqt_intlog2(3*(_i)+1)/2)

static void stack_push(gint *stack, gint *nstack, gint i)

{
  if ( *nstack >= SQT_TRAVERSE_STACK_SIZE )
    g_error("%s: stack overflow, "
	    "increase SQT_TRAVERSE_STACK_SIZE and recompile",
	    __FUNCTION__) ;
  stack[(*nstack)] = i ;
  (*nstack) ++ ;
  
  return ;
}

static void stack_pop(gint *stack, gint *nstack, gint *i)

{
  *i = stack[(*nstack)-1] ;
  (*nstack) -- ;
  
  return ;
}

static void tree_children(gint i, gint c[])

{
  gint d, off, j ;
  /*find depth of i*/
  d = sqt_tree_level_depth(i) ;

  off = (sqt_fourpow(d) - 1)/3 ;
  j = i-off ;

  off = (sqt_fourpow(d+1) - 1)/3 ;
  c[0] = off + 4*j + 0 ;
  c[1] = off + 4*j + 1 ;
  c[2] = off + 4*j + 2 ;
  c[3] = off + 4*j + 3 ;

  return ;
}

static void tree_quad_kw(SQT_REAL *ce, gint ne, gint Nk,
			 SQT_REAL wt,
			 SQT_REAL *q, gint nq, 
#ifdef SQT_SINGLE_PRECISION
			 sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
			 sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
			 SQT_REAL *quad, gint nc,
			 SQT_REAL *work,
			 gpointer data)

{
  gint i, j, ns ;
  SQT_REAL s[4], t[4], w[4], J[4], y[12], n[12], si, ti, *st ;
  gboolean init ;

  st = &(quad[nc]) ;
  ns = 4*(nq/4) ;
  /* init = TRUE ; */
  init = FALSE ;
  for ( i = 0 ; i < ns ; i += 4 ) {
    si = q[3*(i+0)+0] ; ti = q[3*(i+0)+1] ; w[0] = q[3*(i+0)+2] ;
    /*interpolating coordinates on unit triangle*/
    s[0] = (1.0-si-ti)*st[0] + si*st[2] + ti*st[4] ;
    t[0] = (1.0-si-ti)*st[1] + si*st[3] + ti*st[5] ;

    si = q[3*(i+1)+0] ; ti = q[3*(i+1)+1] ; w[1] = q[3*(i+1)+2] ;
    /*interpolating coordinates on unit triangle*/
    s[1] = (1.0-si-ti)*st[0] + si*st[2] + ti*st[4] ;
    t[1] = (1.0-si-ti)*st[1] + si*st[3] + ti*st[5] ;

    si = q[3*(i+2)+0] ; ti = q[3*(i+2)+1] ; w[2] = q[3*(i+2)+2] ;
    /*interpolating coordinates on unit triangle*/
    s[2] = (1.0-si-ti)*st[0] + si*st[2] + ti*st[4] ;
    t[2] = (1.0-si-ti)*st[1] + si*st[3] + ti*st[5] ;

    si = q[3*(i+3)+0] ; ti = q[3*(i+3)+1] ; w[3] = q[3*(i+3)+2] ;
    /*interpolating coordinates on unit triangle*/
    s[3] = (1.0-si-ti)*st[0] + si*st[2] + ti*st[4] ;
    t[3] = (1.0-si-ti)*st[1] + si*st[3] + ti*st[5] ;
    
    /*coordinates on physical triangle*/
    SQT_FUNCTION_NAME(sqt_element_interp_vector)(ce, ne, Nk, s, t, 4,
						 y, 3, n, 3, J, work) ;

    func(s[0], t[0], w[0]*J[0]*wt, &(y[3*0]), &(n[3*0]), &(work[3*ne*0]), ne,
	 quad, nc, init, data) ;
    func(s[1], t[1], w[1]*J[1]*wt, &(y[3*1]), &(n[3*1]), &(work[3*ne*1]), ne,
	 quad, nc, init, data) ;
    func(s[2], t[2], w[2]*J[2]*wt, &(y[3*2]), &(n[3*2]), &(work[3*ne*2]), ne,
	 quad, nc, init, data) ;
    func(s[3], t[3], w[3]*J[3]*wt, &(y[3*3]), &(n[3*3]), &(work[3*ne*3]), ne,
	 quad, nc, init, data) ;
    init = FALSE ;
  }

  ns = nq % 4 ; i = nq - ns ;
  for ( j = 0 ; j < ns ; j ++ ) {
    si = q[3*(i+j)+0] ; ti = q[3*(i+j)+1] ; w[j] = q[3*(i+j)+2] ;
    /*interpolating coordinates on unit triangle*/
    s[j] = (1.0-si-ti)*st[0] + si*st[2] + ti*st[4] ;
    t[j] = (1.0-si-ti)*st[1] + si*st[3] + ti*st[5] ;
    SQT_FUNCTION_NAME(sqt_element_interp)(ce, ne, Nk, s[j], t[j], y, n, J,
  					  NULL, work) ;
    func(s[j], t[j], w[j]*J[0]*wt, y, n, work, ne, quad, nc, 0, data) ;
  }

  return ;
}

static gint tree_traverse(gint idx,
			  SQT_REAL *ce, gint ne, gint Nk,
			  SQT_REAL wt,
			  SQT_REAL *q, gint nq,
#ifdef SQT_SINGLE_PRECISION
			  sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
			  sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
			  SQT_REAL *quad, gint nc,
			  SQT_REAL tol,
			  gpointer data,
			  SQT_REAL *work)

{
  gint d, c[4], w, len, i, nt ;
  SQT_REAL *qp, *qc[4], *st, *stc, err, dq ;

  len = nc + 6 ;
  d = sqt_tree_level_depth(idx) ;
  /*index of current quadrature in quad (parent)*/
  w  = 4*d + (idx - sqt_tree_level_offset(d)) % 4 ;
  qp = &(quad[w*len]) ;  
  tree_children(idx, c) ;
  
  w  = 4*(d+1) + (c[0] - sqt_tree_level_offset(d+1)) % 4 ;
  qc[0] = &(quad[w*len]) ;
  memset(qc[0], 0, 4*len*sizeof(SQT_REAL)) ;
  /*parent element vertices*/
  st = &(qp[nc]) ;
  wt = wt/sqt_fourpow(d+1) ;
  /*split st and integrate on child nodes*/
  for ( i = 0 ; i < 4 ; i ++ ) {
    w  = 4*(d+1) + (c[i] - sqt_tree_level_offset(d+1)) % 4 ;
    qc[i] = &(quad[w*len]) ;
    stc = &((qc[i])[nc]) ; sqt_triangle_divide_loop(i, st, stc) ;
    tree_quad_kw(ce, ne, Nk, wt, q, nq, func, qc[i], nc, work, data) ;
  }

  /*subtract qp and add qc*/
  err = 0.0 ; nt = 0 ;
  for ( i = 0 ; i < nc ; i ++ ) {
    dq = qc[0][i] + qc[1][i] + qc[2][i] + qc[3][i] - qp[i] ;
    quad[i] += dq ;
    err = MAX(fabs(dq), err) ;
    if ( fabs(dq) > tol ) nt ++ ;
  }

  /* fprintf(stderr, "node: %d; depth: %d; nt == %d\n", idx, d, nt) ; */
  
  if ( err > tol ) return 0 ;

  return 1 ;
}

static void tree_subquad_kw(SQT_REAL *ce, gint ne, gint Nk,
			    SQT_REAL *q, gint nq,
#ifdef SQT_SINGLE_PRECISION
			    sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
			    sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
			    SQT_REAL *quad,
			    gint nc,
			    SQT_REAL tol, gint dmax,
			    gpointer data, SQT_REAL *work,
			    SQT_REAL wt, SQT_REAL *st)

{
  SQT_REAL *kwork ;
  gint stack[SQT_TRAVERSE_STACK_SIZE], i, j, s, c[4], nstack ;
  
  /*
   * size of workspace for each integral evaluation:
   *   number of integrands (nc)
   *   triangle vertices (3x2)
   *
   * workspace size: 4*dmax*(nc+6) + 12*ne
   */

  work[nc+0] = st[0] ; work[nc+1] = st[1] ; 
  work[nc+2] = st[2] ; work[nc+3] = st[3] ; 
  work[nc+4] = st[4] ; work[nc+5] = st[5] ;
  kwork = &(work[4*(dmax+2)*(nc+6)]) ;

  tree_quad_kw(ce, ne, Nk, wt, q, nq, func, work, nc, kwork, data) ;
  memset(work, 0 , nc*sizeof(gdouble)) ;
  
  nstack = 0 ;
  stack_push(stack, &nstack, 0) ;
  
  while ( nstack != 0 ) {
    j = stack[nstack-1] ;
    stack_pop(stack, &nstack, &i) ;
    
    s = tree_traverse(j, ce, ne, Nk, wt, q, nq, func, work, nc, tol,
    		      data, kwork) ;

    if ( s == 0 && sqt_tree_level_depth(j) < dmax ) {
      tree_children(j, c) ;
      stack_push(stack, &nstack, c[3]) ;
      stack_push(stack, &nstack, c[2]) ;
      stack_push(stack, &nstack, c[1]) ;
      stack_push(stack, &nstack, c[0]) ;
    }
  }
  
  /* fprintf(stderr, "%lg\n", work[0]) ; */
  /* memcpy(quad, work, nc*sizeof(SQT_REAL)) ; */

  for ( i = 0 ; i < nc ; i  ++ ) quad[i] += work[i] ;
  
  return ;
}

gint SQT_FUNCTION_NAME(sqt_tree_quad_kw)(SQT_REAL *ce, gint ne, gint Nk,
					 SQT_REAL *q, gint nq,
#ifdef SQT_SINGLE_PRECISION
					 sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
					 sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
					 SQT_REAL *quad,
					 gint nc,
					 SQT_REAL tol, gint dmax,
					 gpointer data, SQT_REAL *work)

/*
 * work space size: 4*dmax*nc + 4*ne*3 ;
 */
{
  SQT_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0}, stsub[6], *kwork ;
  gint stack[SQT_TRAVERSE_STACK_SIZE], i, j, s, c[4], nstack ;
  
  /*
   * size of workspace for each integral evaluation:
   *   number of integrands (nc)
   *   triangle vertices (3x2)
   *
   * workspace size: 4*dmax*(nc+6) + 12*ne
   */

  /* tree_subquad_kw(ce, ne, Nk, q, nq, func, quad, nc, */
  /* 		  tol, dmax, data, work, 1.0, st) ; */
  /* return 0 ; */
  
  for ( i = 0 ; i < 4 ; i ++ ) {
    sqt_triangle_divide_loop(i, st, stsub) ;
    tree_subquad_kw(ce, ne, Nk, q, nq, func, quad, nc,
		    tol, dmax, data, work, 0.25, stsub) ;
  }

  /* fprintf(stderr, "%lg\n", quad[0]) ; */

  return 0 ;
  
  work[nc+0] = st[0] ; work[nc+1] = st[1] ; 
  work[nc+2] = st[2] ; work[nc+3] = st[3] ; 
  work[nc+4] = st[4] ; work[nc+5] = st[5] ;
  kwork = &(work[4*(dmax+2)*(nc+6)]) ;

  tree_quad_kw(ce, ne, Nk, 1.0, q, nq, func, work, nc, kwork, data) ;
  memset(work, 0 , nc*sizeof(gdouble)) ;
  
  nstack = 0 ;
  stack_push(stack, &nstack, 0) ;
  
  while ( nstack != 0 ) {
    j = stack[nstack-1] ;
    stack_pop(stack, &nstack, &i) ;
    
    s = tree_traverse(j, ce, ne, Nk, 1.0, q, nq, func, work, nc, tol,
    		      data, kwork) ;

    if ( s == 0 && sqt_tree_level_depth(j) < dmax ) {
      tree_children(j, c) ;
      stack_push(stack, &nstack, c[3]) ;
      stack_push(stack, &nstack, c[2]) ;
      stack_push(stack, &nstack, c[1]) ;
      stack_push(stack, &nstack, c[0]) ;
    }
  }
  
  fprintf(stderr, "%lg\n", work[0]) ;
  memcpy(quad, work, nc*sizeof(SQT_REAL)) ;

  /* fprintf(stderr, "%lg\n", quad[0]) ; */
  
  return 0 ;
}

