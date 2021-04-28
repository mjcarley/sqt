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

#include "config.h"

#include "sqt-private.h"

gint dgesdd_(gchar *JOBZ, gint *M, gint *N, gdouble *A, gint *lda,
	     gdouble *S, gdouble *U, gint *ldu, gdouble *VT, gint *ldvt,
	     gdouble *work, gint *lwork, gint *iwork, gint *info) ;
gint sgesdd_(gchar *JOBZ, gint *M, gint *N, gfloat *A, gint *lda,
	     gfloat *S, gfloat *U, gint *ldu, gfloat *VT, gint *ldvt,
	     gfloat *work, gint *lwork, gint *iwork, gint *info) ;

static SQT_REAL tri_orientation(SQT_REAL *xt)

{
  SQT_REAL det, A[4] ;

  A[0] = xt[1*2+0] - xt[0*2+0] ; A[1] = xt[1*2+1] - xt[0*2+1] ;
  A[2] = xt[2*2+0] - xt[0*2+0] ; A[3] = xt[2*2+1] - xt[0*2+1] ;
  
  det = A[0]*A[3] - A[1]*A[2] ;

  return det ;
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
  SQT_REAL L[10], dLds[10], dLdt[10], n[3], U[9], S[2], V[4] ;
  SQT_REAL work[138], Ai[4], dr[6], tmp ;
  gint i, two = 2, three = 3, lwork, info, iwork ;

  /*physical location of singularity and reference Jacobian*/
  SQT_FUNCTION_NAME(sqt_element_shape_3d)(nt, s, t, L, dLds, dLdt,
					  NULL, NULL, NULL) ;
  SQT_FUNCTION_NAME(sqt_element_point_interp_3d)(xt, xstr, nt, L, dLds, dLdt,
						 x0, n, J0) ;

  memset(dr, 0, 6*sizeof(SQT_REAL)) ;
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
  sgesdd_("A", &three, &two, dr, &three, S, U, &three, V, &two,
	  work, &lwork, &iwork, &info) ;
#else
  /*dr has been filled in FORTRAN ordering so this is a
    standard call to dgesdd with no transposition*/
  dgesdd_("A", &three, &two, dr, &three, S, U, &three, V, &two,
	  work, &lwork, &iwork, &info) ;
  g_assert(info == 0) ;
#endif /*SQT_SINGLE_PRECISION*/
  
  /*A = V*inv(S) = V*[1/S(1) 0; 0 1/S(2)*/
  A[0] = V[0]/S[0] ; A[1] = V[1]/S[1] ;
  A[2] = V[2]/S[0] ; A[3] = V[3]/S[1] ;
  sqt_invert2x2(Ai, A) ;

  /*remapped unit triangle*/
  xti[0] = Ai[0]*(0 - s) + Ai[1]*(0 - t) ;
  xti[1] = Ai[2]*(0 - s) + Ai[3]*(0 - t) ;
  xti[2] = Ai[0]*(1 - s) + Ai[1]*(0 - t) ;
  xti[3] = Ai[2]*(1 - s) + Ai[3]*(0 - t) ;
  xti[4] = Ai[0]*(0 - s) + Ai[1]*(1 - t) ;
  xti[5] = Ai[2]*(0 - s) + Ai[3]*(1 - t) ;

  /*check orientation of remapped triangle*/
  if ( tri_orientation(xti) > 0.0 ) return 0 ;

  tmp = xti[2*0+0] ; xti[2*0+0] = xti[2*1+0] ; xti[2*1+0] = tmp ;
  tmp = xti[2*0+1] ; xti[2*0+1] = xti[2*1+1] ; xti[2*1+1] = tmp ;
  
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
      /* rr = 0.5*d*(1.0 + qr[2*j+0]) ; */
      rr = d*qr[2*j+0] ;
      /*coordinates on unit triangle*/
      s1 = s + rr*Ac ; t1 = t + rr*As ;
      /*map to physical triangle*/
      SQT_FUNCTION_NAME(sqt_element_point_3d)(xt, xstr, nt, s1, t1, y, n, &J) ;
      /*quadrature weight*/
      /* wt = J/J0*rr*qr[2*j+1]*0.5*d*2.0*M_PI/nqt ; */
      /* wt = J/J0*rr*qr[2*j+1]*d*2.0*M_PI/nqt ; */
      wt = J/J0*d*qr[2*j+1]*d*2.0*M_PI/nqt ;
      func(s1, t1, wt, y, n, NULL, 0, Iq, nqi, 0, data) ;
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
	func(s1, t1, wt, y, n, NULL, 0, Iq, nqi, 0, data) ;
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
  SQT_REAL x0[3], J0, xti[6], A[4], d ;

  triangle_mapping(xe, xstr, ne, s0, t0, x0, &J0, A, xti) ;

  d = MIN(point_line_distance2(&(xti[0*2]), &(xti[1*2])),
	  point_line_distance2(&(xti[1*2]), &(xti[2*2]))) ;
  d = MIN(point_line_distance2(&(xti[2*2]), &(xti[0*2])), d) ;

  d = 0.5*SQRT(d) ;
  
  triangle_quad(func, data, xti, xe, xstr, ne, A, s0, t0, J0, d, N,
		quad, nc) ;
  
  return 0 ;
}

static gint triangle_mapping_kw(SQT_REAL *ce, gint ne, gint Nk,
				SQT_REAL s, SQT_REAL t,
				SQT_REAL *x0, SQT_REAL *J0, SQT_REAL *A,
				SQT_REAL *xti, SQT_REAL *work)

/*
  mapping information for transformation of unit triangle to physical
  element via remapped triangle with conformal mapping at singularity

  ce:    interpolation coefficients for nodes of physical triangle
  ne:    number of interpolation coefficients 
  Nk:    order of Koornwinder polynomials
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
  SQT_REAL n[3], U[9], S[2], V[4], Ai[4], dr[6], tmp ;
  gint two = 2, three = 3, lwork, info, iwork ;
  
  /*physical location of singularity and reference Jacobian*/
  SQT_FUNCTION_NAME(sqt_element_interp)(ce, ne, Nk, s, t, x0, n, J0, dr, work) ;
  
  /*lwork = 138 from a workspace query: size of problem is fixed so
    lwork is too*/
  lwork = 138 ; info = 1 ; iwork = 1 ;
#ifdef SQT_SINGLE_PRECISION
  /*need to put single-precision BLAS call in here*/
  g_assert_not_reached() ;
  sgesdd_("A", &three, &two, dr, &three, S, U, &three, V, &two,
	  work, &lwork, &iwork, &info) ;
#else
  /*dr has been filled in FORTRAN ordering so this is a standard call
    to dgesdd with no transposition*/
  dgesdd_("A", &three, &two, dr, &three, S, U, &three, V, &two,
	  work, &lwork, &iwork, &info) ;
  g_assert(info == 0) ;
#endif /*SQT_SINGLE_PRECISION*/
  
  /*A = V*inv(S) = V*[1/S(1) 0; 0 1/S(2)*/
  A[0] = V[0]/S[0] ; A[1] = V[1]/S[1] ;
  A[2] = V[2]/S[0] ; A[3] = V[3]/S[1] ;
  sqt_invert2x2(Ai, A) ;

  /*remapped unit triangle*/
  xti[0] = Ai[0]*(0 - s) + Ai[1]*(0 - t) ;
  xti[1] = Ai[2]*(0 - s) + Ai[3]*(0 - t) ;
  xti[2] = Ai[0]*(1 - s) + Ai[1]*(0 - t) ;
  xti[3] = Ai[2]*(1 - s) + Ai[3]*(0 - t) ;
  xti[4] = Ai[0]*(0 - s) + Ai[1]*(1 - t) ;
  xti[5] = Ai[2]*(0 - s) + Ai[3]*(1 - t) ;

  /*check orientation of remapped triangle*/
  if ( tri_orientation(xti) > 0.0 ) return 0 ;

  tmp = xti[2*0+0] ; xti[2*0+0] = xti[2*1+0] ; xti[2*1+0] = tmp ;
  tmp = xti[2*0+1] ; xti[2*0+1] = xti[2*1+1] ; xti[2*1+1] = tmp ;
  
  return 0 ;
}

static gint triangle_quad_kw(
#ifdef SQT_SINGLE_PRECISION
			  sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
			  sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
			  gpointer data,
			  SQT_REAL *xti,
			  SQT_REAL *ce, gint ne, gint Nk,
			  SQT_REAL *A, SQT_REAL s, SQT_REAL t,
			  SQT_REAL J0,
			  SQT_REAL d, gint N,
			  SQT_REAL *Iq, gint nqi, SQT_REAL *work)

/*
  func: integrand
  data: data for func
  xti:  mapped plane triangle
  ce:   interpolation coefficients for nodes of physical triangle
  ne:   number of interpolation coefficients 
  Nk:   order of Koornwinder polynomials
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
  SQT_REAL Ac, As, Cd, Sd, tmp, c ;
  gint i, j, k, idx[] = {0, 1, 2, 0}, nqr, nqt ;

  /*quadrature on inner disc*/
  nqt = 3*(N+2) + 2 ;
  nqr = (N+1)/2 + 1 ;

  c = d*2.0*M_PI/nqt/J0 ;
  SQT_FUNCTION_NAME(legendre_quadrature_select)(N, &qr, &nqr) ;
  C = 1.0 ; S = 0.0 ;
  Cd = cos(2.0*M_PI/nqt) ; Sd = sin(2.0*M_PI/nqt) ; 
  for ( i = 0 ; i < nqt ; i ++ ) {
    Ac = A[0]*C + A[1]*S ;
    As = A[2]*C + A[3]*S ;
    for ( j = 0 ; j < nqr ; j ++ ) {
      /* rr = 0.5*d*(1.0 + qr[2*j+0]) ; */
      rr = d*qr[2*j+0] ;
      /*coordinates on unit triangle*/
      s1 = s + rr*Ac ; t1 = t + rr*As ;
      /*map to physical triangle*/
      SQT_FUNCTION_NAME(sqt_element_interp)(ce, ne, Nk, s1, t1, y, n, &J,
					    NULL, work) ;
      /*quadrature weight*/
      /* wt = J*rr*qr[2*j+1]*c ; */
      wt = J*d*qr[2*j+1]*c ;
      func(s1, t1, wt, y, n, work, ne, Iq, nqi, 0, data) ;
    }
    tmp = C ;
    C = C*Cd - S*Sd ;
    S = tmp*Sd + S*Cd ;
  }

  for ( k = 0 ; k < 3 ; k ++ ) {
    SQT_REAL S0 ;
    /*singularity is at origin on remapped triangle*/
    triangle_decomp(zero, &(xti[2*idx[k]]), &(xti[2*idx[k+1]]),
		    &r0, &th0, &r, &th) ;
    sgn = SIGN(th0) ; th0 = fabs(th0) ;
    S0 = sin(th0) ;

    c = th0/J0 ;
    /*select th quadrature*/
    SQT_FUNCTION_NAME(angular_quadrature_select)(N, r0, th0, &qt, &nqt) ;
    for ( i = 0 ; i < nqt ; i ++ ) {
      ti = th0*qt[i*2+0] ;
      C = cos(sgn*ti+th) ; S = sin(sgn*ti+th) ;
      Ac = A[0]*C + A[1]*S ;
      As = A[2]*C + A[3]*S ;
      /*radial limit and quadrature selection*/
      M = r*r0*S0/(r0*sin(th0-ti) + sin(ti)) ;
      SQT_FUNCTION_NAME(radial_quadrature_select)(N, d/M, &qr, &nqr) ;
      for ( j = 0 ; j < nqr ; j ++ ) {
      	rr = d + (M-d)*qr[2*j+0] ;
	/*coordinates on unit triangle*/
	s1 = s + rr*Ac ; t1 = t + rr*As ;
	/*map to physical triangle*/
	SQT_FUNCTION_NAME(sqt_element_interp)(ce, ne, Nk, s1, t1, y, n, &J,
					      NULL, work) ;
	/*quadrature weight*/
	wt = (M-d)*qr[2*j+1]*qt[2*i+1]*rr*J*c ;
	func(s1, t1, wt, y, n, work, ne, Iq, nqi, 0, data) ;
      }
    }
  }
  
  return 0 ;
}

gint SQT_FUNCTION_NAME(sqt_singular_quad_kw)(SQT_REAL *ce, gint ne, gint Nk,
					     SQT_REAL s0, SQT_REAL t0,
					     gint N,
#ifdef SQT_SINGLE_PRECISION
					     sqt_quadrature_func_f_t func,
#else /*SQT_SINGLE_PRECISION*/
					     sqt_quadrature_func_t func,
#endif /*SQT_SINGLE_PRECISION*/
					     SQT_REAL *quad, gint nc,
					     gpointer data,
					     SQT_REAL *work)
/*workspace size 12*ne*/
  
{
  SQT_REAL x0[3], J0, xti[64], A[4], d ;

  triangle_mapping_kw(ce, ne, Nk, s0, t0, x0, &J0, A, xti, work) ;

  d = MIN(point_line_distance2(&(xti[0*2]), &(xti[1*2])),
	  point_line_distance2(&(xti[1*2]), &(xti[2*2]))) ;
  d = MIN(point_line_distance2(&(xti[2*2]), &(xti[0*2])), d) ;

  d = 0.5*SQRT(d) ;
  
  triangle_quad_kw(func, data, xti, ce, ne, Nk, A, s0, t0, J0, d, N,
		   quad, nc, work) ;
  
  return 0 ;
}
