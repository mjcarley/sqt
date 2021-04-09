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

/*basic integration routine used for testing and well-separated elements*/

gint SQT_FUNCTION_NAME(sqt_basic_quad_tri)(SQT_REAL *xe, gint xstr, gint ne,
					   SQT_REAL *q, gint nq,
#ifdef SQT_SINGLE_PRECISION
					   sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
					   sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
					   SQT_REAL *quad, gint nc,
					   gpointer data)


{
  SQT_REAL s, t, w, J, y[3], n[3] ;
  gint i ;
  
  memset(quad, 0, nc*sizeof(SQT_REAL)) ;

  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[3*i+0] ; t = q[3*i+1] ; w = q[3*i+2] ;
    SQT_FUNCTION_NAME(sqt_element_point_3d)(xe, xstr, ne, s, t, y, n, &J) ;
    w *= J ;
    func(s, t, w, y, n, NULL, 0, quad, nc, 0, data) ;
  }

  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_basic_quad_kw)(SQT_REAL *xe, gint xstr,
					  SQT_REAL *n , gint nstr,
					  SQT_REAL *wt, gint wstr,
					  gint ne,
#ifdef SQT_SINGLE_PRECISION
					   sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
					   sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
					   SQT_REAL *quad, gint nc,
					   gpointer data)


{
  gint i ;
  
  memset(quad, 0, nc*sizeof(SQT_REAL)) ;

  g_assert_not_reached() ;
  for ( i = 0 ; i < ne ; i ++ ) {
    func(0, 0, wt[i*wstr], &(xe[i*xstr]), &(n [i*nstr]), NULL, 0,
	 quad, nc, 0, data) ;
  }

  return 0 ;
}
