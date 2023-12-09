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

inline static void shape3(SQT_REAL s, SQT_REAL t, SQT_REAL L[])

{
  L[0] = 1.0 - s - t ; 
  L[1] =       s     ; 
  L[2] =           t ;

  return ;
}

static gint calc_point(SQT_REAL *xe, gint xstr, gint ne,
		       SQT_REAL *L, SQT_REAL *p)

{
  gint i ;

  p[0] = p[1] = p[2] = 0.0 ;

  for ( i = 0 ; i < ne ; i ++ ) {
    p[0] += L[i]*xe[i*xstr+0] ;
    p[1] += L[i]*xe[i*xstr+1] ;
    p[2] += L[i]*xe[i*xstr+2] ;
  }
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_element_shape_3d)(gint ne, SQT_REAL s, SQT_REAL t,
					     SQT_REAL *L,
					     SQT_REAL *dLds, SQT_REAL *dLdt,
					     SQT_REAL *dLdss, SQT_REAL *dLdst,
					     SQT_REAL *dLdtt)

{
  switch ( ne ) {
  default: g_assert_not_reached() ; break ;
  case 3:
    L[0] = 1.0 - s - t ; 
    L[1] =       s     ; 
    L[2] =           t ;
    if ( dLds == NULL ) break ;
    dLds[0] = -1.0 ; dLds[1] =  1.0 ; dLds[2] =  0.0 ;
    dLdt[0] = -1.0 ; dLdt[1] =  0.0 ; dLdt[2] =  1.0 ;
    if ( dLdss == NULL ) break ;
    dLdss[0] = 0.0 ; dLdss[1] =  0.0 ; dLdss[2] =  0.0 ;
    dLdst[0] = 0.0 ; dLdst[1] =  0.0 ; dLdst[2] =  0.0 ;
    dLdtt[0] = 0.0 ; dLdtt[1] =  0.0 ; dLdtt[2] =  0.0 ;
    break ;
  case 6:
    L[0] = 2.0*(1.0-s-t)*(0.5-s-t) ;
    L[1] = 2.0*s*(s-0.5) ;
    L[2] = 2.0*t*(t-0.5) ;
    L[3] = 4.0*s*(1.0-s-t) ;
    L[4] = 4.0*s*t ;
    L[5] = 4.0*t*(1.0-s-t) ;

    if ( dLds == NULL ) break ;

    dLds[0] = -3.0 + 4.0*s + 4.0*t ;
    dLds[1] = -1.0 + 4.0*s ;
    dLds[2] =  0.0 ;
    dLds[3] =  4.0 - 8.0*s - 4.0*t ;
    dLds[4] =  4.0*t ;
    dLds[5] = -4.0*t ;

    dLdt[0] = -3.0 + 4.0*s + 4.0*t ;
    dLdt[1] =  0.0 ;
    dLdt[2] = -1.0 + 4.0*t ;
    dLdt[3] = -4.0*s ;
    dLdt[4] =  4.0*s ;
    dLdt[5] =  4.0 - 4.0*s - 8.0*t ;

    if ( dLdss == NULL ) break ;

    dLdss[0] =  4.0 ;
    dLdss[1] =  4.0 ;
    dLdss[2] =  0.0 ;
    dLdss[3] = -8.0 ;
    dLdss[4] =  0.0 ;
    dLdss[5] =  0.0 ;

    dLdst[0] =  4.0 ;
    dLdst[1] =  0.0 ;
    dLdst[2] =  0.0 ;
    dLdst[3] = -4.0 ;
    dLdst[4] =  4.0 ;
    dLdst[5] = -4.0 ;

    dLdtt[0] =  4.0 ;
    dLdtt[1] =  0.0 ;
    dLdtt[2] =  4.0 ;
    dLdtt[3] =  0.0 ;
    dLdtt[4] =  0.0 ;
    dLdtt[5] = -8.0 ;

    break ;
  }

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_element_point_interp_3d)(SQT_REAL *xe,
						    gint xstr, gint ne,
						    SQT_REAL *L,
						    SQT_REAL *dLds,
						    SQT_REAL *dLdt,
						    SQT_REAL *y,
						    SQT_REAL *n,
						    SQT_REAL *J)

{
  gint i ;
  SQT_REAL dyds[3]={0.0}, dydt[3]={0.0} ;
  
  y[0] = y[1] = y[2] = n[0] = n[1] = n[2] = 0.0 ;
  for ( i = 0 ; i < ne ; i ++ ) {
    y[0] += xe[i*xstr+0]*L[i] ; y[1] += xe[i*xstr+1]*L[i] ;
    y[2] += xe[i*xstr+2]*L[i] ; 
    dyds[0] += xe[i*xstr+0]*dLds[i] ; dyds[1] += xe[i*xstr+1]*dLds[i] ;
    dyds[2] += xe[i*xstr+2]*dLds[i] ; 
    dydt[0] += xe[i*xstr+0]*dLdt[i] ; dydt[1] += xe[i*xstr+1]*dLdt[i] ;
    dydt[2] += xe[i*xstr+2]*dLdt[i] ; 
  }

  sqt_vector_cross(n, dyds, dydt) ;
  
  *J = sqt_vector_length(n) ;

  n[0] /= (*J) ; n[1] /= (*J) ; n[2] /= (*J) ;

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_element_point_3d)(SQT_REAL *xe, gint xstr, gint ne,
					     SQT_REAL s, SQT_REAL t,
					     SQT_REAL *y, SQT_REAL *n,
					     SQT_REAL *J)

{
  SQT_REAL L[32]={0.0}, dLds[32]={0.0}, dLdt[32]={0.0} ;
  
  SQT_FUNCTION_NAME(sqt_element_shape_3d)(ne, s, t, L, dLds, dLdt,
					  NULL, NULL, NULL) ;
  SQT_FUNCTION_NAME(sqt_element_point_interp_3d)(xe, xstr, ne, L,
						 dLds, dLdt, y, n, J) ;
  
  return 0 ;
}

SQT_REAL SQT_FUNCTION_NAME(sqt_element_area)(SQT_REAL *xe, gint xstr, gint ne,
					     SQT_REAL *qrule, gint nq)

{
  gint i ;
  SQT_REAL y[3], s, t, w, n[3], J, A ;

  A = 0.0 ;
  for ( i = 0 ; i < nq ; i ++ ) {
    s = qrule[3*i+0] ; t = qrule[3*i+1] ; w = qrule[3*i+2] ; 
    SQT_FUNCTION_NAME(sqt_element_point_3d)(xe, xstr, ne, s, t, y, n, &J) ;
    A += w*J ;
  }
  
  return A ;
}

gint SQT_FUNCTION_NAME(sqt_quadrature_select)(gint nq, SQT_REAL **q,
					      gint *order)

{
#ifdef SQT_SINGLE_PRECISION
  switch ( nq ) {
  default:
    g_error("%s: %d point quadrature not available "
	    "(available rules are 7, 25, 54, 85, 126, 175, 453)",
	    __FUNCTION__, nq) ;
  case   7: *q = WANDZURA_7_F   ; *order =  5 ; break ;
  case  25: *q = WANDZURA_25_F  ; *order = 10 ; break ;
  case  54: *q = WANDZURA_54_F  ; *order = 15 ; break ;
  case  85: *q = WANDZURA_85_F  ; *order = 20 ; break ;
  case 126: *q = WANDZURA_126_F ; *order = 25 ; break ;
  case 175: *q = WANDZURA_175_F ; *order = 30 ; break ;
  case 453: *q = XIAO_GIMBUTAS_453_F ; *order = 50 ; break ;
  }
#else /*SQT_SINGLE_PRECISION*/
  switch ( nq ) {
  default:
    g_error("%s: %d point quadrature not available "
	    "(available rules are 7, 25, 54, 85, 126, 175, 453)",
	    __FUNCTION__, nq) ;
    break ;
  case   7: *q = WANDZURA_7   ; *order =  5 ; break ;
  case  25: *q = WANDZURA_25  ; *order = 10 ; break ;
  case  54: *q = WANDZURA_54  ; *order = 15 ; break ;
  case  85: *q = WANDZURA_85  ; *order = 20 ; break ;
  case 126: *q = WANDZURA_126 ; *order = 25 ; break ;
  case 175: *q = WANDZURA_175 ; *order = 30 ; break ;
  case 453: *q = XIAO_GIMBUTAS_453 ; *order = 50 ; break ;
  }
#endif /*SQT_SINGLE_PRECISION*/
  
  return 0 ;
}

/* gint SQT_FUNCTION_NAME(sqt_quadrature_optimal_points)(SQT_REAL Rb, */
/* 						      SQT_REAL r0, */
/* 						      SQT_REAL r1, */
/* 						      gint nr, */
/* 						      gint order, gint ngp, */
/* 						      gint pmax, gint smax, */
/* 						      SQT_REAL tol, */
/* 						      SQT_REAL *x0, */
/* 						      SQT_REAL *n, */
/* 						      SQT_REAL *x, */
/* 						      SQT_REAL *rc, */
/* 						      gint *pq, gint *s) */

/* { */
/*   gint ncalc, ncmin, pc, sc, i ; */
/*   SQT_REAL r, rp, c[3] ; */
  
/*   ncmin = G_MAXINT ; */

/*   for ( i = 0 ; i <= nr ; i ++ ) { */
/*     rp = r0 + (r1-r0)*i/nr ; */
/*     sqt_vector_shift(c,x0,n,rp) ; */
/*     r = sqt_vector_distance(x,c) ; */
/*     SQT_FUNCTION_NAME(sqt_truncation_optimal)(Rb, r, rp, order, pmax, smax, */
/* 					      tol, &pc, &sc) ; */
/*     ncalc = ngp*(1 << sc) ; */
/*     if ( ncalc < ncmin ) { */
/*       ncmin = ncalc ; */
/*       *pq = pc ; *s = sc ; */
/*       *rc = rp ; */
/*     } */
/*   } */

/*   return 0 ; */
/* } */

/* gint SQT_FUNCTION_NAME(sqt_triangle_expansion)(SQT_REAL *xe, gint xstr, gint ne, */
/* 					       SQT_REAL *fe, gint fstr, gint nf, */
/* 					       SQT_REAL *q, gint nq, gint oq, */
/* 					       SQT_REAL *x, SQT_REAL tol, */
/* 					       gint Nmax, gint smax, */
/* 					       SQT_REAL *c, SQT_REAL *C, */
/* 					       gint *N, gint *s) */

/* { */
/*   /\* SQT_REAL Ae, rb ; *\/ */

/*   /\* Ae = SQT_FUNCTION_NAME(sqt_element_area)(xe, xstr, ne, q, nq) ; *\/ */
/*   /\* rb = SQRT(4.0*Ae/M_PI) ; *\/ */

/*   return 0 ; */
/* } */
						  
gint SQT_FUNCTION_NAME(sqt_triangle_curvature)(SQT_REAL *xe, gint xstr, gint ne,
					       SQT_REAL s, SQT_REAL t,
					       SQT_REAL *kg, SQT_REAL *km)

{
  SQT_REAL p1[3], ps[3], pt[3] ;
  SQT_REAL L0[16], Ls[16], Lt[16], Lss[16], Lst[16], Ltt[16],
    ds, n[3], pss[3], pst[3], ptt[3] ;
  SQT_REAL E, F, G, L, M, N ;

  SQT_FUNCTION_NAME(sqt_element_shape_3d)(ne, s, t, L0, Ls, Lt,
					  Lss, Lst, Ltt) ;
  calc_point(xe, xstr, ne, L0,  p1) ;
  calc_point(xe, xstr, ne, Ls, ps) ;
  calc_point(xe, xstr, ne, Lt, pt) ;
  calc_point(xe, xstr, ne, Lss, pss) ;
  calc_point(xe, xstr, ne, Lst, pst) ;
  calc_point(xe, xstr, ne, Ltt, ptt) ;

  sqt_vector_cross(n,ps,pt) ;
  ds = sqt_vector_length(n) ;
  n[0] /= ds ; n[1] /= ds ; n[2] /= ds ;
  
  /*first fundamental forms*/
  E = sqt_vector_scalar(ps, ps) ;
  F = sqt_vector_scalar(ps, pt) ;
  G = sqt_vector_scalar(pt, pt) ;

  /*second fundamental forms*/
  L = sqt_vector_scalar(pss, n) ;
  M = sqt_vector_scalar(pst, n) ;
  N = sqt_vector_scalar(ptt, n) ;

  *kg = (L*N - M*M)/(E*G - F*F) ;
  *km = 0.5*(E*N - 2*F*M + G*L)/(E*G - F*F) ;
  
  return 0 ;
}
