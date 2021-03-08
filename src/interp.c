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

  xs = &(buf[3]) ; xt = &(buf[6]) ;

#ifndef SQT_SINGLE_PRECISION
  gint i3 = 3 ;
  SQT_REAL al = 1.0, bt = 0.0 ;

  K = work ; Ks = &(K[nq]) ; Kt = &(Ks[nq]) ;

  sqt_koornwinder_deriv_nm(Nk, s, t, K, 1, Ks, 1, Kt, 1, nq) ;
  blaswrap_dgemm(FALSE, FALSE, i3, i3, nq, al, K, nq, ci, i3, bt, buf, i3) ;

  x[0] = buf[0] ; x[1] = buf[1] ; x[2] = buf[2] ; 
  
  sqt_vector_cross(n, xs, xt) ;
  *J = sqt_vector_length(n) ;
  n[0] /= (*J) ; n[1] /= (*J) ; n[2] /= (*J) ; 
#else
  g_assert_not_reached() ;
#endif /*SQT_SINGLE_PRECISION*/

  if ( dx == NULL ) return 0 ;

  memcpy(dx, &(buf[3]), 6*sizeof(SQT_REAL)) ;
  
  /* dx[0] = xs[0] ; dx[1] = xs[1] ; dx[2] = xs[2] ; */
  /* dx[3] = xt[0] ; dx[4] = xt[1] ; dx[5] = xt[2] ; */
  
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
#endif /*SQT_SINGLE_PRECISION*/
  }
  
  return 0 ;
}
