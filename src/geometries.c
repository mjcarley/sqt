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

#include "config.h"

#include "sqt-private.h"

gint SQT_FUNCTION_NAME(sqt_geometry_stellarator)(SQT_REAL u, SQT_REAL v,
						 SQT_REAL *x, SQT_REAL *n)

/*
  stellarator geometry from Greengard et al.
*/
  
{
  SQT_REAL d[] = {-1, -1,  0.17,
		 -1,  0,  0.11,
		 +0,  0,  1.00,
		 1, 0, 4.5,
		 2, 0, -0.25,
		 0, 1, 0.07,
		 2, 1, -0.45} ;
  SQT_REAL i, j, dij, xu[3], xv[3] ;
  gint k ;

  x [0] = x [1] = x [2] = 0.0 ;
  xu[0] = xu[1] = xu[2] = 0.0 ;
  xv[0] = xv[1] = xv[2] = 0.0 ;

  for ( k = 0 ; k < 7 ; k ++ ) {
    i = d[3*k+0] ; j = d[3*k+1] ; dij = d[3*k+2] ;
    x[0] += dij*cos(v)*cos((1.0-i)*u + j*v) ;
    x[1] += dij*sin(v)*cos((1.0-i)*u + j*v) ;
    x[2] += dij*       sin((1.0-i)*u + j*v) ;

    xu[0] += -(1.0-i)*dij*cos(v)*sin((1.0-i)*u + j*v) ;
    xu[1] += -(1.0-i)*dij*sin(v)*sin((1.0-i)*u + j*v) ;
    xu[2] +=  (1.0-i)*dij*       cos((1.0-i)*u + j*v) ;

    xv[0] +=
      dij*(-sin(v)*cos((1.0-i)*u + j*v) - j*cos(v)*sin((1.0-i)*u + j*v)) ;
    xv[1] +=
      dij*(cos(v)*cos((1.0-i)*u + j*v) - j*sin(v)*sin((1.0-i)*u + j*v)) ;
    xv[2] +=
      dij*j*cos((1.0-i)*u + j*v) ;
  }

  sqt_vector_cross(n, xu, xv) ;
  i = sqt_vector_length(n) ;
  n[0] /= i ; n[1] /= i ; n[2] /= i ;
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_geometry_sphere)(SQT_REAL th, SQT_REAL ph,
					    SQT_REAL *x, SQT_REAL *xu,
					    SQT_REAL *xv)

{
  x[0] = cos(th)*sin(ph) ;
  x[1] = sin(th)*sin(ph) ;
  x[2] =         cos(ph) ;

  xu[0] = -sin(th)*sin(ph) ;
  xu[1] =  cos(th)*sin(ph) ;
  xu[2] =  0.0 ;

  xv[0] = cos(th)*cos(ph) ;
  xv[1] = sin(th)*cos(ph) ;
  xv[2] =        -sin(ph) ;
    
  return 0 ;
}

