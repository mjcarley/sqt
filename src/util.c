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
#include <math.h>

#include <glib.h>

#include <sqt.h>

#include "sqt-private.h"

gint SQT_FUNCTION_NAME(sqt_cartesian_to_spherical)(SQT_REAL *x0,
						   SQT_REAL *x,
						   SQT_REAL *r,
						   SQT_REAL *th,
						   SQT_REAL *ph)

{
  *r = sqt_vector_distance2(x, x0) ;
  if ( *r == 0.0 ) { *ph = *th = 0.0 ; return 0 ; }

  *r = SQRT((*r)) ;
  *ph = ATAN2(x[1]-x0[1], x[0]-x0[0]) ;

  *th = ACOS((x[2]-x0[2])/(*r)) ;

  return 0 ;
}
