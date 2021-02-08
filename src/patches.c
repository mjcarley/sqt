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

gint SQT_FUNCTION_NAME(sqt_patch_nodes_tri)(SQT_REAL *xe, gint xstr, gint ne,
					    SQT_REAL *s, gint sstr,
					    SQT_REAL *t, gint tstr,
					    gint nst,
					    SQT_REAL *xp, gint pstr,
					    SQT_REAL *np, gint nstr)

{
  gint i ;
  SQT_REAL n[3], J ;
  
  g_assert(ne == 3 || ne == 6) ;

  if ( np == NULL ) {
    np = n ; nstr = 0 ;
  }
  
  for ( i = 0 ; i < nst ; i ++ ) {
    sqt_element_point_3d(xe, xstr, ne, s[i*sstr], t[i*tstr],
			 &(xp[i*pstr]), &(np[i*nstr]), &J) ;  
  }
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_patch_nodes_sphere)(SQT_REAL th0, SQT_REAL ph0,
					       SQT_REAL th1, SQT_REAL ph1,
					       SQT_REAL th2, SQT_REAL ph2,
					       SQT_REAL *s, gint sstr,
					       SQT_REAL *t, gint tstr,
					       gint nst,
					       SQT_REAL *xp, gint pstr,
					       SQT_REAL *np, gint nstr)

{
  gint i ;
  SQT_REAL n[3], J, th, ph, xv[3], xu[3] ;
  
  if ( np == NULL ) {
    np = n ; nstr = 0 ;
  }
  
  for ( i = 0 ; i < nst ; i ++ ) {
    th = th0*(1.0 - s[i*sstr] - t[i*tstr]) + th1*s[i*sstr] + th2*t[i*tstr] ;
    ph = ph0*(1.0 - s[i*sstr] - t[i*tstr]) + ph1*s[i*sstr] + ph2*t[i*tstr] ;

    sqt_geometry_sphere(th, ph, &(xp[i*pstr]), xu, xv) ;
    np[i*nstr+0] = xp[i*pstr+0] ;
    np[i*nstr+1] = xp[i*pstr+1] ;
    np[i*nstr+2] = xp[i*pstr+2] ;
  }
  
  return 0 ;
}
