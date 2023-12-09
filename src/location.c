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

gint SQT_FUNCTION_NAME(sqt_element_nearest_point)(SQT_REAL *xe, gint xstr,
						  gint ne,
						  SQT_REAL *xf, 
						  SQT_REAL *sn, SQT_REAL *tn, 
						  SQT_REAL *p1, SQT_REAL tol,
						  gint nimax, SQT_REAL *kn)

/*
  implementation of Li X, Wu Z, Pan F et al. A geometric strategy
  algorithm for orthogonal projection onto a parametric
  surface. Journal of computer science and technology 34(6): 1279-1293
  Nov. 2019. https://dx.doi.org/10.1007/s11390-019-1967-z
*/
  
{
  SQT_REAL ps[3], pt[3], L0[16], Ls[16], Lt[16], Lss[16], Lst[16], Ltt[16],
    ds, dt, n[3], pss[3], pst[3], ptt[3] ;
  SQT_REAL E, F, G, L, M, N, C4, C5 ;
  gint ni ;

  ni = 0 ;

  do {
    SQT_FUNCTION_NAME(sqt_element_shape_3d)(ne, *sn, *tn, L0, Ls, Lt,
					    NULL, NULL, NULL) ;

    calc_point(xe, xstr, ne, L0, p1) ;
    calc_point(xe, xstr, ne, Ls, ps) ;
    calc_point(xe, xstr, ne, Lt, pt) ;

    /* vector_cross(n, ps, pt) ; */
    E = sqt_vector_scalar(ps, ps) ;
    F = sqt_vector_scalar(ps, pt) ;
    G = sqt_vector_scalar(pt, pt) ;

    C4 = sqt_vector_diff_scalar(xf, p1, ps) ;
    C5 = sqt_vector_diff_scalar(xf, p1, pt) ;

    dt = E*G - F*F ;
    ds =  (C4*G - C5*F)/dt ;
    dt = -(C4*F - C5*E)/dt ;

    *sn += ds ; *tn += dt ;
    
    ni ++ ;
  } while ( ds*ds + dt*dt > tol*tol && ni < nimax ) ;
  
  SQT_FUNCTION_NAME(sqt_element_shape_3d)(ne, *sn, *tn, L0, Ls, Lt,
					  Lss, Lst, Ltt) ;
  calc_point(xe, xstr, ne, L0,  p1) ;
  calc_point(xe, xstr, ne, Ls, ps) ;
  calc_point(xe, xstr, ne, Lt, pt) ;
  calc_point(xe, xstr, ne, Lss, pss) ;
  calc_point(xe, xstr, ne, Lst, pst) ;
  calc_point(xe, xstr, ne, Ltt, ptt) ;

  sqt_vector_cross(n,ps,pt) ;
  ds = sqt_vector_length(n) ;
    /* SQRT(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]) ; */
  n[0] /= ds ; n[1] /= ds ; n[2] /= ds ;
  
  /*first fundamental forms*/
  E = sqt_vector_scalar(ps, ps) ;
  F = sqt_vector_scalar(ps, pt) ;
  G = sqt_vector_scalar(pt, pt) ;

  /*second fundamental forms*/
  L = sqt_vector_scalar(pss, n) ;
  M = sqt_vector_scalar(pst, n) ;
  N = sqt_vector_scalar(ptt, n) ;

  *kn = (L*N - M*M)/(E*G - F*F) ;
  
  return ni ;
}
