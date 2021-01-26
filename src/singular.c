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

#include <sqt.h>

#include "sqt-private.h"

gint dgesdd_(gchar *JOBZ, gint *M, gint *N, gdouble *A, gint *lda,
	     gdouble *S, gdouble *U, gint *ldu, gdouble *VT, gint *ldvt,
	     gdouble *work, gint *lwork, gint *iwork, gint *info) ;


static gint invert2x2(SQT_REAL *A, SQT_REAL *Ai)

{
  SQT_REAL det ;

  det = A[0]*A[3] - A[1]*A[2] ;

  Ai[0] =  A[3]/det ; Ai[1] = -A[1]/det ;
  Ai[2] = -A[2]/det ; Ai[3] =  A[0]/det ;

  return 0 ;
}


static gint triangle_mapping(SQT_REAL *xt, gint xstr, gint nt,
			     SQT_REAL s, SQT_REAL t,
			     SQT_REAL *x0, SQT_REAL *J0, SQT_REAL *A,
			     SQT_REAL *xti)


/*
  mapping information for transformation of unit triangle to physical
  element via remapped triangle with conformal mapping at singularity

  xt:    nodes of physical triangle
  xstr:  stride in xt
  nt:    number of nodes on xt (3 or 6)
  s, t:  coordinates of singular point on unit triangle
  x0:    calculated physical location of singular point
  J0:    reference Jacobian at singularity
  A:     transformation matrix for Affine mapping
  xti:   triangle for integration using singular quadrature rules

  integrate on xti and map xti to coordinates (s1,t1) on unit triangle
  with

  (s1,t1) = [A]*r[\cos\theta \sin\theta]' + (s,t)'

  map to physical triangle using shape functions as usual

  x = L(s1,t1)*[xt]
*/
  
{
  SQT_REAL L[16], dLds[16], dLdt[16], n[16], U[32], S[32], V[32] ;
  SQT_REAL work[138], Ai[16], dr[16] ;
  gint i, one = 1, two = 2, three = 3, lwork, info, iwork ;

  /* memset(L, 0, 16*sizeof(SQT_REAL)) ; */
  /* memset(dLds, 0, 16*sizeof(SQT_REAL)) ; */
  /* memset(dLdt, 0, 16*sizeof(SQT_REAL)) ; */
  /* memset(U, 0, 16*sizeof(SQT_REAL)) ; */
  /* memset(S, 0, 16*sizeof(SQT_REAL)) ; */
  /* memset(V, 0, 16*sizeof(SQT_REAL)) ; */
  memset(dr, 0, 16*sizeof(SQT_REAL)) ;
  memset(Ai, 0, 4*sizeof(SQT_REAL)) ;
  memset(work, 0, 138*sizeof(SQT_REAL)) ;
  memset(n, 0, 3*sizeof(SQT_REAL)) ;

  /*physical location of singularity and reference Jacobian*/
  SQT_FUNCTION_NAME(sqt_element_shape_3d)(nt, s, t, L, dLds, dLdt,
					  NULL, NULL, NULL) ;
  SQT_FUNCTION_NAME(sqt_element_point_interp_3d)(xt, xstr, nt, L, dLds, dLdt,
						 x0, n, J0) ;

  /* dr[2] = dr[3] = dr[4] = dr[5] = 0.0 ; */
  
  /*SVD of Jacobian matrix*/
  for ( i = 0 ; i < nt ; i ++ ) {
    dr[0] += dLds[i]*xt[i*xstr+0] ;
    dr[1] += dLds[i]*xt[i*xstr+1] ;
    dr[2] += dLds[i]*xt[i*xstr+2] ;
    dr[3] += dLdt[i]*xt[i*xstr+0] ;
    dr[4] += dLdt[i]*xt[i*xstr+1] ;
    dr[5] += dLdt[i]*xt[i*xstr+2] ;
  }

  /*lwork = 138 from a workspace query: size of problem is fixed so
    lwork is too*/
  lwork = 138 ; info = 1 ; iwork = 1 ;
#ifdef SQT_SINGLE_PRECISION
  /*need to put single-precision BLAS call in here*/
  g_assert_not_reached() ;
#else
  dgesdd_("A", &three, &two, dr, &three, S, U, &three, V, &two,
	  work, &lwork, &iwork, &info) ;
  g_assert(info == 0) ;
#endif /*SQT_SINGLE_PRECISION*/
  
  /*A = V*inv(S) = V*[1/S(1) 0; 0 1/S(2)*/
  A[0] = V[0]/S[0] ; A[1] = V[1]/S[1] ;
  A[2] = V[2]/S[0] ; A[3] = V[3]/S[1] ;
  invert2x2(A, Ai) ;

  /*remapped unit triangle*/
  xti[0] = Ai[0]*(0 - s) + Ai[1]*(0 - t) ;
  xti[1] = Ai[2]*(0 - s) + Ai[3]*(0 - t) ;
  xti[2] = Ai[0]*(1 - s) + Ai[1]*(0 - t) ;
  xti[3] = Ai[2]*(1 - s) + Ai[3]*(0 - t) ;
  xti[4] = Ai[0]*(0 - s) + Ai[1]*(1 - t) ;
  xti[5] = Ai[2]*(0 - s) + Ai[3]*(1 - t) ;
  
  return 0 ;
}

static gdouble atan2_uw(gdouble y, gdouble x)

{
  gdouble th ;

  th = atan2(y, x) ;

  if ( th >= 0 ) return th ;

  return 2.0*M_PI + th ;
}

static gint triangle_decomp(SQT_REAL *x0, SQT_REAL *x1, SQT_REAL *x2,
			    SQT_REAL *r0, SQT_REAL *th0, SQT_REAL *r,
			    SQT_REAL *th)

{
  SQT_REAL r1, r2 ;

  r1 = sqrt((x0[0] - x1[0])*(x0[0] - x1[0]) +
	    (x0[1] - x1[1])*(x0[1] - x1[1])) ;
  r2 = sqrt((x0[0] - x2[0])*(x0[0] - x2[0]) +
	    (x0[1] - x2[1])*(x0[1] - x2[1])) ;

  if ( r1 > r2 ) {
    *r = r1 ; *r0 = r2/r1 ;
    *th = atan2_uw(x1[1] - x0[1], x1[0] - x0[0]) ;
    *th0 = atan2_uw(x2[1] - x0[1], x2[0] - x0[0]) - *th ;

    if ( *th0 < -M_PI ) *th0 += 2.0*M_PI ;

    return 0 ;
  }

  *r = r2 ; *r0 = r1/r2 ;
  *th = atan2_uw(x2[1] - x0[1], x2[0] - x0[0]) ;
  *th0 = atan2_uw(x1[1] - x0[1], x1[0] - x0[0]) - *th ;

  if ( *th0 > 0 ) *th0 -= 2.0*M_PI ;
  
  return 0 ;
}

static gint triangle_quad(
#ifdef SQT_SINGLE_PRECISION
			  sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
			  sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
			  gpointer data,
			  SQT_REAL *xti,
			  SQT_REAL *xt, gint xstr, gint nt,
			  SQT_REAL *A, SQT_REAL s, SQT_REAL t,
			  SQT_REAL J0,
			  SQT_REAL d, gint N,
			  SQT_REAL *Iq, gint nqi)

/*
  func: integrand
  data: data for func
  xti:  mapped plane triangle
  xt:   physical space triangle
  xstr: stride in xt
  nt:   number of nodes on xt (3 or 6)
  A:    mapping matrix
  s,t:  singularity coordinates on unit triangle
  J0:   Jacobian at singularity
  d:    radius of central disc
  N:    order of integration
  Iq:   integrals of func over triangle
  nqi:  number of functions to integrate
*/
  
{
  SQT_REAL r0, th0, r, th, *qr, *qt, rr, ti, M, sgn ;
  SQT_REAL y[3], s1, t1, J, n[3], wt, C, S, zero[]={0,0} ;
  SQT_REAL Ac, As ;
  gint i, j, k, idx[] = {0, 1, 2, 0}, nqr, nqt ;

  /* memset(Iq, 0, nqi*sizeof(SQT_REAL)) ; */
  
  /*quadrature on inner disc*/
  nqt = 3*(N+2) + 2 ;
  nqr = (N+1)/2 + 1 ;

  SQT_FUNCTION_NAME(legendre_quadrature_select)(N, &qr, &nqr) ;
  for ( i = 0 ; i < nqt ; i ++ ) {
    ti = 2.0*M_PI*i/nqt ;
    C = cos(ti) ; S = sin(ti) ;
    Ac = A[0]*C + A[1]*S ;
    As = A[2]*C + A[3]*S ;
    for ( j = 0 ; j < nqr ; j ++ ) {
      rr = 0.5*d*(1.0 + qr[2*j+0]) ;
      /*coordinates on unit triangle*/
      s1 = s + rr*Ac ; t1 = t + rr*As ;
      /*map to physical triangle*/
      SQT_FUNCTION_NAME(sqt_element_point_3d)(xt, xstr, nt, s1, t1, y, n, &J) ;
      /*quadrature weight*/
      wt = J/J0*rr*qr[2*j+1]*0.5*d*2.0*M_PI/nqt ;
      func(s1, t1, wt, y, n, Iq, nqi, data) ;
    }
  }

  for ( k = 0 ; k < 3 ; k ++ ) {
    /*singularity is at origin on remapped triangle*/
    triangle_decomp(zero, &(xti[2*idx[k]]), &(xti[2*idx[k+1]]),
		    &r0, &th0, &r, &th) ;
    sgn = SIGN(th0) ; th0 = fabs(th0) ;
    /*select th quadrature*/
    SQT_FUNCTION_NAME(angular_quadrature_select)(N, r0, th0, &qt, &nqt) ;
    for ( i = 0 ; i < nqt ; i ++ ) {
      ti = th0*qt[i*2+0] ;
      C = cos(sgn*ti+th) ; S = sin(sgn*ti+th) ;
      Ac = A[0]*C + A[1]*S ;
      As = A[2]*C + A[3]*S ;
      /*radial limit and quadrature selection*/
      M = r0*sin(th0)/(r0*sin(th0-ti) + sin(ti)) ;
      SQT_FUNCTION_NAME(radial_quadrature_select)(N, d/M/r, &qr, &nqr) ;
      for ( j = 0 ; j < nqr ; j ++ ) {
      	rr = d + (M*r-d)*qr[2*j+0] ;
	/*coordinates on unit triangle*/
	s1 = s + rr*Ac ; t1 = t + rr*As ;
	/*map to physical triangle*/
	SQT_FUNCTION_NAME(sqt_element_point_3d)(xt, xstr, nt, s1, t1,
						y, n, &J) ;
	/*quadrature weight*/
	wt = (M*r-d)*qr[2*j+1]*qt[2*i+1]*rr*th0*J/J0 ;
	func(s1, t1, wt, y, n, Iq, nqi, data) ;
      }
    }
  }
  
  return 0 ;
}

static SQT_REAL point_line_distance2(SQT_REAL *x0, SQT_REAL *x1)

{
  SQT_REAL l2, r[2], x[2], t ;

  r[0] = x1[0] - x0[0] ; r[1] = x1[1] - x0[1] ;
  l2 = r[0]*r[0] + r[1]*r[1] ;

  t = -(x0[0]*r[0] + x0[1]*r[1])/l2 ;
  x[0] = x0[0] + t*r[0] ;
  x[1] = x0[1] + t*r[1] ;
  
  l2 = x[0]*x[0] + x[1]*x[1] ;

  return l2 ;
}

gint SQT_FUNCTION_NAME(sqt_singular_quad_tri)(SQT_REAL *xe, gint xstr, gint ne,
					      SQT_REAL s0, SQT_REAL t0,
					      gint N,
#ifdef SQT_SINGLE_PRECISION
					      sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
					      sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
					      SQT_REAL *quad, gint nc,
					      gpointer data)
  
{
  SQT_REAL x0[3], J0, xti[64], A[4], d ;

  triangle_mapping(xe, xstr, ne, s0, t0, x0, &J0, A, xti) ;

  d = MIN(point_line_distance2(&(xti[0*2]), &(xti[1*2])),
	  point_line_distance2(&(xti[1*2]), &(xti[2*2]))) ;
  d = MIN(point_line_distance2(&(xti[2*2]), &(xti[0*2])), d) ;

  d = 0.5*SQRT(d) ;
  
  triangle_quad(func, data, xti, xe, xstr, ne, A, s0, t0, J0, d, N,
		quad, nc) ;
  
  return 0 ;
}
