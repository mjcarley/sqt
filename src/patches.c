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

gint SQT_FUNCTION_NAME(sqt_patch_nodes_tri)(SQT_REAL *xe, gint xstr, gint ne,
					    SQT_REAL *s, gint sstr,
					    SQT_REAL *t, gint tstr,
					    SQT_REAL *w, gint wstr,
					    gint nst,
					    SQT_REAL *xp, gint pstr,
					    SQT_REAL *np, gint nstr,
					    SQT_REAL *wt, gint wtstr)

{
  gint i ;
  SQT_REAL n[3], J, buffer[] = {1.0, 1.0} ;
  
  g_assert(ne == 3 || ne == 6) ;

  if ( np == NULL ) {
    np = n ; nstr = 0 ;
  }

  if ( w == NULL ) {
    w = &(buffer[0]) ; wstr = 0 ;
  }

  if ( wt == NULL ) {
    wt = &(buffer[1]) ; wtstr = 0 ;
  }
  
  for ( i = 0 ; i < nst ; i ++ ) {
    SQT_FUNCTION_NAME(sqt_element_point_3d)(xe, xstr, ne, s[i*sstr], t[i*tstr],
					    &(xp[i*pstr]), &(np[i*nstr]), &J) ;
    wt[i*wtstr] = J*w[i*wstr] ;
  }
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_patch_nodes_sphere)(SQT_REAL rho,
					       SQT_REAL th0, SQT_REAL ph0,
					       SQT_REAL th1, SQT_REAL ph1,
					       SQT_REAL th2, SQT_REAL ph2,
					       SQT_REAL *s, gint sstr,
					       SQT_REAL *t, gint tstr,
					       SQT_REAL *w, gint wstr,
					       gint nst,
					       SQT_REAL *xp, gint pstr,
					       SQT_REAL *np, gint nstr,
					       SQT_REAL *wt, gint wtstr)

{
  gint i ;
  SQT_REAL n[3], J, th, ph, xph[3], xth[3], buffer[] = {1.0, 1.0} ;
  SQT_REAL xs[3], xt[3], ths, tht, phs, pht ;
    
  if ( np == NULL ) {
    np = n ; nstr = 0 ;
  }
  
  if ( w == NULL ) {
    w = &(buffer[0]) ; wstr = 0 ;
  }

  if ( wt == NULL ) {
    wt = &(buffer[1]) ; wtstr = 0 ;
  }

  for ( i = 0 ; i < nst ; i ++ ) {
    th = th0*(1.0 - s[i*sstr] - t[i*tstr]) + th1*s[i*sstr] + th2*t[i*tstr] ;
    ph = ph0*(1.0 - s[i*sstr] - t[i*tstr]) + ph1*s[i*sstr] + ph2*t[i*tstr] ;

    ths = th0*(0.0 - 1.0 - 0.0) + th1*1.0 + th2*0.0 ;
    tht = th0*(0.0 - 0.0 - 1.0) + th1*0.0 + th2*1.0 ;
    phs = ph0*(0.0 - 1.0 - 0.0) + ph1*1.0 + ph2*0.0 ;
    pht = ph0*(0.0 - 0.0 - 1.0) + ph1*0.0 + ph2*1.0 ;
    
    SQT_FUNCTION_NAME(sqt_geometry_sphere)(th, ph, &(xp[i*pstr]), xth, xph) ;
    xp[i*pstr+0] *= rho ; xp[i*pstr+1] *= rho ; xp[i*pstr+2] *= rho ; 

    np[i*nstr+0] = xp[i*pstr+0] ;
    np[i*nstr+1] = xp[i*pstr+1] ;
    np[i*nstr+2] = xp[i*pstr+2] ;

    xs[0] = xth[0]*ths + xph[0]*phs ;
    xs[1] = xth[1]*ths + xph[1]*phs ;
    xs[2] = xth[2]*ths + xph[2]*phs ;
    xt[0] = xth[0]*tht + xph[0]*pht ;
    xt[1] = xth[1]*tht + xph[1]*pht ;
    xt[2] = xth[2]*tht + xph[2]*pht ;

    sqt_vector_cross(n, xs, xt) ;

    J = sqt_vector_length(n)*rho*rho ;
    
    wt[i*wtstr] = J*w[i*wstr] ;
  }
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_patch_nodes_ellipsoid)(SQT_REAL a, SQT_REAL b,
						  SQT_REAL c,
						  SQT_REAL th0, SQT_REAL ph0,
						  SQT_REAL th1, SQT_REAL ph1,
						  SQT_REAL th2, SQT_REAL ph2,
						  SQT_REAL *s, gint sstr,
						  SQT_REAL *t, gint tstr,
						  SQT_REAL *w, gint wstr,
						  gint nst,
						  SQT_REAL *xp, gint pstr,
						  SQT_REAL *np, gint nstr,
						  SQT_REAL *wt, gint wtstr)

{
  gint i ;
  SQT_REAL n[3], J, th, ph, xph[3], xth[3], buffer[] = {1.0, 1.0} ;
  SQT_REAL xs[3], xt[3], ths, tht, phs, pht ;
    
  if ( np == NULL ) {
    np = n ; nstr = 0 ;
  }
  
  if ( w == NULL ) {
    w = &(buffer[0]) ; wstr = 0 ;
  }

  if ( wt == NULL ) {
    wt = &(buffer[1]) ; wtstr = 0 ;
  }

  for ( i = 0 ; i < nst ; i ++ ) {
    th = th0*(1.0 - s[i*sstr] - t[i*tstr]) + th1*s[i*sstr] + th2*t[i*tstr] ;
    ph = ph0*(1.0 - s[i*sstr] - t[i*tstr]) + ph1*s[i*sstr] + ph2*t[i*tstr] ;

    ths = th0*(0.0 - 1.0 - 0.0) + th1*1.0 + th2*0.0 ;
    tht = th0*(0.0 - 0.0 - 1.0) + th1*0.0 + th2*1.0 ;
    phs = ph0*(0.0 - 1.0 - 0.0) + ph1*1.0 + ph2*0.0 ;
    pht = ph0*(0.0 - 0.0 - 1.0) + ph1*0.0 + ph2*1.0 ;
    
    SQT_FUNCTION_NAME(sqt_geometry_sphere)(th, ph, &(xp[i*pstr]), xth, xph) ;
    xp[i*pstr+0] *= a ; xp[i*pstr+1] *= b ; xp[i*pstr+2] *= c ; 
    xth[0] *= a ; xth[1] *= b ; xth[2] *= c ; 
    xph[0] *= a ; xph[1] *= b ; xph[2] *= c ; 

    xs[0] = xth[0]*ths + xph[0]*phs ;
    xs[1] = xth[1]*ths + xph[1]*phs ;
    xs[2] = xth[2]*ths + xph[2]*phs ;
    xt[0] = xth[0]*tht + xph[0]*pht ;
    xt[1] = xth[1]*tht + xph[1]*pht ;
    xt[2] = xth[2]*tht + xph[2]*pht ;

    sqt_vector_cross(n, xs, xt) ;

    J = sqt_vector_length(n) ;

    np[i*nstr+0] = n[0]/J ;
    np[i*nstr+1] = n[1]/J ;
    np[i*nstr+2] = n[2]/J ;
    
    wt[i*wtstr] = J*w[i*wstr] ;
  }
  
  return 0 ;
}
