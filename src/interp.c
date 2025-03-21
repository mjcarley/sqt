/* This file is part of SQT, a library for Singular Quadrature on Triangles
 *
 * Copyright (C) 2021 Michael Carley
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

#include <blaswrap.h>

#include <sqt.h>

#include "sqt-private.h"

gint SQT_FUNCTION_NAME(sqt_element_interp)(SQT_REAL *ci, gint nq,
					   gint Nk,
					   SQT_REAL s, SQT_REAL t,
					   SQT_REAL *x, SQT_REAL *n,
					   SQT_REAL *J, SQT_REAL *dx,
					   SQT_REAL *work)

  /*
   * ci is found by multiplying element nodes by Koornwinder matrix:
   * al = 1.0 ; bt = 0.0 ;
   * blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr, bt, ci, i3) ;
   * 
   * outputs x: interpolated position
   *         n: surface normal at x
   *         J: Jacobian determinant
   *        dx: if not NULL, Jacobian matrix
   *
   * workspace must be at least 3*nq elements for calculation of
   * Koornwinder polynomials and derivatives
   */

{
  SQT_REAL *K, *Ks, *Kt, buf[9], *xs, *xt ;

  g_assert(work != NULL) ;
  
  xs = &(buf[3]) ; xt = &(buf[6]) ;

  gint i3 = 3 ;
  SQT_REAL al = 1.0, bt = 0.0 ;

  K = work ; Ks = &(K[nq]) ; Kt = &(Ks[nq]) ;

  SQT_FUNCTION_NAME(sqt_koornwinder_deriv_nm)(Nk, s, t, K, 1, Ks, 1, Kt, 1,
					      nq) ;
#ifndef SQT_SINGLE_PRECISION
  blaswrap_dgemm(FALSE, FALSE, i3, i3, nq, al, K, nq, ci, i3, bt, buf, i3) ;
#else
  g_assert_not_reached() ;
  blaswrap_sgemm(FALSE, FALSE, i3, i3, nq, al, K, nq, ci, i3, bt, buf, i3) ;
#endif /*SQT_SINGLE_PRECISION*/

  x[0] = buf[0] ; x[1] = buf[1] ; x[2] = buf[2] ; 
  
  sqt_vector_cross(n, xs, xt) ;
  *J = sqt_vector_length(n) ;
  n[0] /= (*J) ; n[1] /= (*J) ; n[2] /= (*J) ; 

  if ( dx == NULL ) return 0 ;

  memcpy(dx, &(buf[3]), 6*sizeof(SQT_REAL)) ;
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_element_interp_vector)(SQT_REAL *ci, gint nq,
						  gint Nk,
						  SQT_REAL *s, SQT_REAL *t,
						  gint nst,
						  SQT_REAL *x, gint xstr,
						  SQT_REAL *n, gint nstr,
						  SQT_REAL *J,
						  SQT_REAL *work)

  /*
   * ci is found by multiplying element nodes by Koornwinder matrix:
   * al = 1.0 ; bt = 0.0 ;
   * blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr, bt, ci, i3) ;
   * 
   * outputs x: interpolated position
   *         n: surface normal at x
   *         J: Jacobian determinant
   *        dx: if not NULL, Jacobian matrix
   *
   * workspace must be at least 12*nq elements for calculation of
   * Koornwinder polynomials and derivatives for four elements in
   * vectorized calculation
   */

{
  SQT_REAL *K, *Ks, *Kt, buf[36] ;

  gint i3 = 3, nr, i ;
  SQT_REAL al = 1.0, bt = 0.0 ;

  K = work ; Ks = &(K[nq]) ; Kt = &(Ks[nq]) ;
  nr = 3*nst ;
  
  SQT_FUNCTION_NAME(sqt_koornwinder_deriv_nm_vector)(Nk, s, t, nst,
						     K , 1, 3*nq,
						     Ks, 1, 3*nq,
						     Kt, 1, 3*nq, nq) ;
#ifndef SQT_SINGLE_PRECISION
  blaswrap_dgemm(FALSE, FALSE, nr, i3, nq, al, K, nq, ci, i3, bt, buf, i3) ;
#else
  g_assert_not_reached() ;
  blaswrap_sgemm(FALSE, FALSE, nr, i3, nq, al, K, nq, ci, i3, bt, buf, i3) ;
#endif /*SQT_SINGLE_PRECISION*/

  for ( i = 0 ; i < nst ; i ++ ) {
    x[i*xstr+0] = buf[9*i+0] ;
    x[i*xstr+1] = buf[9*i+1] ;
    x[i*xstr+2] = buf[9*i+2] ;
    sqt_vector_cross(&(n[i*nstr]), &(buf[9*i+3]), &(buf[9*i+6])) ;
    J[i] = sqt_vector_length(&(n[i*nstr])) ;
    n[i*nstr+0] /= J[i] ;
    n[i*nstr+1] /= J[i] ;
    n[i*nstr+2] /= J[i] ;
  }
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_interp_matrix)(SQT_REAL *K, gint nk, gint Nk,
					  SQT_REAL *s, gint sstr,
					  SQT_REAL *t, gint tstr,
					  gint nst,
					  SQT_REAL *Ki, SQT_REAL *work)

{
  SQT_REAL al, bt, *Knm ;
  gint i, i1 = 1 ;
  
  Knm = work ;

  al = 1.0 ; bt = 0.0 ;
  for ( i = 0 ; i < nst ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_koornwinder_nm)(Nk, s[i*sstr], t[i*tstr],
					  Knm, 1, nk) ;
#ifndef SQT_SINGLE_PRECISION
    blaswrap_dgemv(TRUE, nk, nk, al, K, nk, Knm, i1, bt, &(Ki[i*nk]), i1) ;
#else /*SQT_SINGLE_PRECISION*/
    g_assert_not_reached() ;
    blaswrap_sgemv(TRUE, nk, nk, al, K, nk, Knm, i1, bt, &(Ki[i*nk]), i1) ;
#endif /*SQT_SINGLE_PRECISION*/
  }
  
  return 0 ;
}
