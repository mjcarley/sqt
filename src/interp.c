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
   * workspace must be at least 3*(Nk+1)*(Nk+2)/2 elements for calculation
   * of Koornwinder polynomials and derivatives
   */

{
  SQT_REAL *K, *Ks, *Kt, xs[3], xt[3] ;
  gint str = 3 ;
  
  K = work ; Ks = &(work[1]) ; Kt = &(work[2]) ;

  sqt_koornwinder_deriv_nm(Nk, s, t, K, str, Ks, str, Kt, str) ;
#ifndef SQT_SINGLE_PRECISION
  gint i1 = 1, i3 = 3 ;
  SQT_REAL al = 1.0, bt = 0.0 ;
  
  blaswrap_dgemv(TRUE, nq, i3, al, ci, i3, K , str, bt, x , i1) ;
  blaswrap_dgemv(TRUE, nq, i3, al, ci, i3, Ks, str, bt, xs, i1) ;
  blaswrap_dgemv(TRUE, nq, i3, al, ci, i3, Kt, str, bt, xt, i1) ;

  sqt_vector_cross(n, xs, xt) ;
  *J = sqt_vector_length(n) ;
  n[0] /= (*J) ; n[1] /= (*J) ; n[2] /= (*J) ; 
#else
  g_assert_not_reached() ;
#endif /*SQT_SINGLE_PRECISION*/

  if ( dx == NULL ) return 0 ;

  dx[0] = xs[0] ; dx[1] = xs[1] ; dx[2] = xs[2] ;
  dx[3] = xt[0] ; dx[4] = xt[1] ; dx[5] = xt[2] ;
  
  return 0 ;
}
