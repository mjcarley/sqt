/* This file is part of SQT, a library for Quadrature By Expansion
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
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <glib.h>

#include <sqt.h>

#include <blaswrap.h>

#include "sqt-private.h"

gchar *tests[] = {"closest_point",
		  "koornwinder",
		  "koornwinder_orthogonality",
		  "blas",
		  "koornwinder_interpolation",
		  "adaptive",
		  "normal",
		  "weights",
		  "matrix_adaptive",
		  "matrix_self",
		  "element_interpolation",
		  "spherical",
		  "koornwinder_derivatives",
		  ""} ;

GTimer *timer ;

static gint invert2x2(SQT_REAL *A, SQT_REAL *Ai)

{
  SQT_REAL det ;

  det = A[0]*A[3] - A[1]*A[2] ;

  Ai[0] =  A[3]/det ; Ai[1] = -A[1]/det ;
  Ai[2] = -A[2]/det ; Ai[3] =  A[0]/det ;

  det = A[0]*A[3] - A[1]*A[2] ;

  return 0 ;
}

static gint parse_test(gchar *arg)

{
  gint i = 0 ;
  
  while ( strlen(tests[i]) != 0) {
    if ( !strcmp(tests[i], arg) ) return i ;
    i ++ ;
  }

  return -1 ;
}

static gint read_element(FILE *f, gdouble *xe, gint *ne, gint *xstr)

{
  gint i ;

  *xstr = 4 ;
  fscanf(f, "%d", ne) ;

  for ( i = 0 ; i < (*ne) ; i ++ ) {
    fscanf(f, "%lg", &(xe[i*(*xstr)+0])) ;
    fscanf(f, "%lg", &(xe[i*(*xstr)+1])) ;
    fscanf(f, "%lg", &(xe[i*(*xstr)+2])) ;
  }
  
  return 0 ;
}


static gint element_closest_point_test(gdouble *xe, gint xstr, gint ne,
				       gdouble s0,  gdouble t0,
				       gdouble *x0, gdouble rc)

{
  gdouble xc[3], n[3], xf[3], J, s, t, xn[3], r[3], kn, km ;
  gint ni ;
  
  fprintf(stderr, "closest point test\n") ;
  fprintf(stderr, "===================\n") ;

  fprintf(stderr, "element: %d nodes\n", ne) ;
  fprintf(stderr, "rc: %lg\n", rc) ;
  fprintf(stderr, "s,t: %lg %lg\n", s0, t0) ;

  /*closest point on element for test purposes*/
  sqt_element_point_3d(xe, xstr, ne, s0, t0, xc, n, &J) ;  
  fprintf(stderr, "xc: %lg %lg %lg\n", xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "n:  %lg %lg %lg\n", n[0], n[1], n[2]) ;
  sqt_triangle_curvature(xe, xstr, ne, s0, t0, &kn, &km) ;
  fprintf(stderr, "kn: %lg\n", kn) ;
  
  xf[0] = xc[0] + n[0]*rc ;
  xf[1] = xc[1] + n[1]*rc ;
  xf[2] = xc[2] + n[2]*rc ;

  fprintf(stderr, "xf: %lg %lg %lg\n", xf[0], xf[1], xf[2]) ;

  /* s = 1/3.0 ; t = 1/3.0 ; */
  /* s = t = 0.0 ; */
  s = 1.0 ; t = 0.0 ;
  /* s = s0 ; t = t0+0.000001 ; */
  ni = sqt_element_nearest_point(xe, xstr, ne, xf, &s, &t, xn, 1e-9, 256, &kn) ;
  
  sqt_element_point_3d(xe, xstr, ne, s, t, xc, n, &J) ;  

  r[0] = xf[0] - xc[0] ; r[1] = xf[1] - xc[1] ; r[2] = xf[2] - xc[2] ;
  fprintf(stderr, "s,t: %lg %lg\n", s, t) ;
  fprintf(stderr, "ni:  %d\n", ni) ;
  fprintf(stderr, "x0:  %lg %lg %lg\n", xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "kn:  %lg\n", kn) ;
  fprintf(stderr, "n:   %lg %lg %lg\n", n[0], n[1], n[2]) ;
  fprintf(stderr, "r:   %lg %lg %lg\n", r[0], r[1], r[2]) ;
  fprintf(stderr, "r.n: %lg (%lg)\n",
	  sqt_vector_scalar(r,n), sqt_vector_length(r)) ;
  
  return 0 ;
}
				       
static gint koornwinder_test(gint N, gdouble u, gdouble v)

{
  gdouble Knm[1024] ;
  gint n, m, str ;

  fprintf(stderr, "koornwinder test\n") ;
  fprintf(stderr, "================\n") ;

  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "(u,v) = (%lg,%lg)\n", u, v) ;

  str = 3 ;
  
  sqt_koornwinder_nm(N, u, v, Knm, str) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    for ( m = 0 ; m <= n ; m ++ ) {
      fprintf(stdout,
	      "%d %d %1.16e %1.16e %1.16e\n",
	      n, m, u, v, Knm[str*(n*(n+1)/2+m)]) ;
    }
  }

  return 0 ;
}

static gint koornwinder_derivatives_test(gint N, gdouble u, gdouble v)

{
  gdouble K[1024], Ku[1024], Kv[1024], eK, eu, ev, ee ;
  gdouble K0[1024], Kup1[1024], Kum1[1024], Kvp1[1024], Kvm1[1024] ;
  gint n, m, str, kstr, ustr, vstr, idx ;

  fprintf(stderr, "koornwinder derivatives test\n") ;
  fprintf(stderr, "============================\n") ;

  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "(u,v) = (%lg,%lg)\n", u, v) ;

  ee = 1e-6 ;
  
  str = 3 ;

  kstr = 4 ; ustr = 2 ; vstr = 3 ;
  
  sqt_koornwinder_nm(N, u, v, K0, str) ;

  sqt_koornwinder_nm(N, u+0.5*ee, v       , Kup1, str) ;
  sqt_koornwinder_nm(N, u-0.5*ee, v       , Kum1, str) ;
  sqt_koornwinder_nm(N, u       , v+0.5*ee, Kvp1, str) ;
  sqt_koornwinder_nm(N, u       , v-0.5*ee, Kvm1, str) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    for ( m = 0 ; m <= n ; m ++ ) {
      idx = n*(n+1)/2 + m ;
      Kup1[str*idx] = (Kup1[str*idx] - Kum1[str*idx])/ee ;
      Kvp1[str*idx] = (Kvp1[str*idx] - Kvm1[str*idx])/ee ;
    }
  }
  
  sqt_koornwinder_deriv_nm(N, u, v, K, kstr, Ku, ustr, Kv, vstr) ;

  eK = eu = ev = 0.0 ;
  for ( n = 0 ; n <= N ; n ++ ) {
    for ( m = 0 ; m <= n ; m ++ ) {
      idx = n*(n+1)/2 + m ;
      eK = MAX(eK, fabs(K0[str*idx]-K[kstr*idx])) ;
      eu = MAX(eu, fabs(Kup1[str*idx]-Ku[ustr*idx])) ;
      ev = MAX(ev, fabs(Kvp1[str*idx]-Kv[vstr*idx])) ;
      fprintf(stdout,
      	      "%d %d %lg %lg %lg %lg %lg %lg %lg %lg\n",
      	      n, m, u, v,
	      K0[str*idx], K[kstr*idx],
	      Kup1[str*idx], Ku[ustr*idx],
	      Kvp1[str*idx], Kv[vstr*idx]) ;
    }
  }

  fprintf(stderr, "error K : %lg\n", eK) ;
  fprintf(stderr, "error Ku: %lg\n", eu) ;
  fprintf(stderr, "error Kv: %lg\n", ev) ;
  
  return 0 ;
}

static gint koornwinder_orthogonality_test(gint N)

{
  gint nq, order, i, n1, m1, n2, m2, idx1, idx2, str ;
  gdouble s, t, Knm[4096], *q, I, w, tol ;

  str = 3 ;
  
  tol = 1e-12 ;
  
  fprintf(stderr, "koornwinder orthogonality test\n") ;
  fprintf(stderr, "==============================\n") ;

  fprintf(stderr, "N = %d\n", N) ;

  nq = 453 ;
  
  sqt_quadrature_select(nq, &q, &order) ;

  for ( n1 = 0 ; n1 <= N ; n1 ++ ) {
    for ( m1 = 0 ; m1 <= n1 ; m1 ++ ) {
      for ( n2 = 0 ; n2 <= N ; n2 ++ ) {
	for ( m2 = 0 ; m2 <= n2 ; m2 ++ ) {
	  idx1 = n1*(n1+1)/2 + m1 ; 
	  idx2 = n2*(n2+1)/2 + m2 ; 
  
	  I = 0.0 ;
	  for ( i = 0 ; i < nq ; i ++ ) {
	    s = q[3*i+0] ; 
	    t = q[3*i+1] ;
	    w = q[3*i+2] ; 
	    
	    sqt_koornwinder_nm(N, s, t, Knm, str) ;
	    
	    I += Knm[str*idx1]*Knm[str*idx2]*w ;
	  }

	  fprintf(stderr, "%d %d %d %d %e ", n1, m1, n2, m2, I) ;
	  if ( n1 != n2 || m1 != m2 ) {
	    if ( fabs(I) > tol ) 
	      fprintf(stderr, "FAIL\n") ;
	    else
	      fprintf(stderr, "PASS\n") ;
	  } else {
	    if ( fabs(I-1.0) > tol ) 
	      fprintf(stderr, "FAIL\n") ;
	    else
	      fprintf(stderr, "PASS\n") ;
	  }	  
	}
      }
    }
  }
  return 0 ;
}

static gint koornwinder_interpolation_test(gint N, gint nq)

{
  gint order, i, i1 = 1 ;
  gdouble s, t, Knm[32768], K[4*65536], *q, f, fr, fi[512], al, bt, c[512] ;
  
  fprintf(stderr, "koornwinder interpolation test\n") ;
  fprintf(stderr, "==============================\n") ;
  fprintf(stderr, "nq = %d\n", nq) ;
  
  sqt_quadrature_select(nq, &q, &order) ;

  N = sqt_koornwinder_interp_matrix(q, nq, K) ;

  fprintf(stderr, "Knm N max: %d\n", N) ;
  
  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[i*3+0] ; t = q[i*3+1] ;
    /* fi[i] = 3.0*s*t - t*t ; */
    fi[i] = sin(2.0*M_PI*s*t/8) ;
  }

  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemv(FALSE, nq, nq, al, K, nq, fi, i1, bt, c, i1) ;

  for ( s = 0 ; s <= 1.0 ; s += 0.1 ) {
    for ( t = 0 ; t <= 1.0-s ; t += 0.1 ) {
      sqt_koornwinder_nm(N, s, t, Knm, 1) ;

      f = blaswrap_ddot(nq, c, i1, Knm, i1) ;
      /* fr = 3.0*s*t - t*t ; */
      fr = sin(2.0*M_PI*s*t/8) ;
  
      fprintf(stderr, "%lg %lg %lg %lg (%lg)\n", s, t, fr, f, fabs(fr-f)) ;
    }
  }

  return 0 ;
}

static gint stellarator(gdouble u, gdouble v, gdouble *x,
			gdouble *xu, gdouble *xv)

/*
  stellarator geometry from Greengard et al.
*/
  
{
  gdouble d[] = {-1, -1,  0.17,
		 -1,  0,  0.11,
		 +0,  0,  1.00,
		 1, 0, 4.5,
		 2, 0, -0.25,
		 0, 1, 0.07,
		 2, 1, -0.45} ;
  gdouble i, j, dij ;
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
  
  return 0 ;
}

static gint element_interpolation_test(gint N, gint nq)

{
  gint order, i, i1 = 1, i2 = 2, i3 = 3, xstr ;
  gdouble s, t, Knm[32768], K[4*65536], *q, f, fr, fi[512], al, bt, c[512] ;
  gdouble Ks[32768], Kt[32768] ;
  gdouble ui[3], vi[3], xi[453*4], ci[453*3], u, v, x[3], xu[3], xv[3] ;
  gdouble y[3], yu[3], yv[3], ys[3], yt[3], us, ut, A[4], Ai[4] ;
  gdouble emax, eumax, evmax ;
  
  fprintf(stderr, "element interpolation test\n") ;
  fprintf(stderr, "==========================\n") ;
  fprintf(stderr, "nq = %d\n", nq) ;
  
  sqt_quadrature_select(nq, &q, &order) ;

  N = sqt_koornwinder_interp_matrix(q, nq, K) ;

  fprintf(stderr, "Knm N max: %d\n", N) ;

  ui[0] = 0.3 ; ui[1] = 0.4  ; ui[2] = 0.35 ;
  vi[0] = 0.1 ; vi[1] = 0.1 ; vi[2] = 0.2 ;

  xstr = 4 ;
  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[i*3+0] ; t = q[i*3+1] ;

    u = ui[0]*(1.0 - s - t) + ui[1]*s + ui[2]*t ;
    v = vi[0]*(1.0 - s - t) + vi[1]*s + vi[2]*t ;

    stellarator(u, v, &(xi[i*xstr]), xu, xv) ;
  }

  /*compute the interpolation coefficients*/
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr, bt, ci, i3) ;

  emax = eumax = evmax = 0.0 ;
  for ( s = 0.01 ; s <= 1.0 ; s += 0.05 ) {
    for ( t = 0.01 ; t <= 1.0 - s ; t += 0.05 ) {
  /* { s = 0.3 ; { t = 0.1 ; */
      u = ui[0]*(1.0 - s - t) + ui[1]*s + ui[2]*t ;
      v = vi[0]*(1.0 - s - t) + vi[1]*s + vi[2]*t ;

      A[0] = ui[1] - ui[0] ; A[1] = vi[1] - vi[0] ;
      A[2] = ui[2] - ui[0] ; A[3] = vi[2] - vi[0] ;
      
      invert2x2(A, Ai) ;
      
      stellarator(u, v, x, xu, xv) ;
      /* fprintf(stdout, "%lg %lg %lg ", x[0], x[1], x[2]) ; */
      /* fprintf(stdout, "%lg %lg %lg ", xu[0], xu[1], xu[2]) ; */
      fprintf(stdout, "%lg %lg %lg ", xv[0], xv[1], xv[2]) ;

      /*interpolate using K*/
      /* sqt_koornwinder_nm(N, s, t, Knm, 1) ; */
      sqt_koornwinder_deriv_nm(N, s, t, Knm, 1, Ks, 1, Kt, 1) ;
      blaswrap_dgemv(TRUE, nq, i3, al, ci, i3, Knm, i1, bt, y , i1) ;
      blaswrap_dgemv(TRUE, nq, i3, al, ci, i3, Ks , i1, bt, ys, i1) ;
      blaswrap_dgemv(TRUE, nq, i3, al, ci, i3, Kt , i1, bt, yt, i1) ;

      yu[0] = Ai[0]*ys[0] + Ai[1]*yt[0] ;
      yu[1] = Ai[0]*ys[1] + Ai[1]*yt[1] ;
      yu[2] = Ai[0]*ys[2] + Ai[1]*yt[2] ;
      
      yv[0] = Ai[2]*ys[0] + Ai[3]*yt[0] ;
      yv[1] = Ai[2]*ys[1] + Ai[3]*yt[1] ;
      yv[2] = Ai[2]*ys[2] + Ai[3]*yt[2] ;
      
      /* fprintf(stdout, "%lg %lg %lg\n", y[0], y[1], y[2]) ; */
      /* fprintf(stdout, "%lg %lg %lg\n", yu[0], yu[1], yu[2]) ; */
      fprintf(stdout, "%lg %lg %lg\n", yv[0], yv[1], yv[2]) ;
      emax = MAX(emax,
		 (x[0]-y[0])*(x[0]-y[0]) +
		 (x[1]-y[1])*(x[1]-y[1]) +
		 (x[2]-y[2])*(x[2]-y[2])) ;
      eumax = MAX(eumax,
		  (xu[0]-yu[0])*(xu[0]-yu[0]) +
		  (xu[1]-yu[1])*(xu[1]-yu[1]) +
		  (xu[2]-yu[2])*(xu[2]-yu[2])) ;
      evmax = MAX(evmax,
		  (xv[0]-yv[0])*(xv[0]-yv[0]) +
		  (xv[1]-yv[1])*(xv[1]-yv[1]) +
		  (xv[2]-yv[2])*(xv[2]-yv[2])) ;
    }
  }

  fprintf(stderr, "maximum interpolation error: %lg\n", sqrt(emax)) ;
  fprintf(stderr, "maximum differentiation error: %lg\n", sqrt(eumax)) ;
  fprintf(stderr, "maximum differentiation error: %lg\n", sqrt(evmax)) ;
  
  return 0 ;
}

static gint sphere(gdouble u, gdouble v, gdouble *x)

{
  x[0] = cos(u)*sin(v) ;
  x[1] = sin(u)*sin(v) ;
  x[2] =        cos(v) ;
  
  return 0 ;
}

static gint spherical_patch_test(gint N, gint nq)

{
  gint order, i, i1 = 1, i2 = 2, i3 = 3, xstr ;
  gdouble s, t, Knm[32768], K[4*65536], *q, f, fr, fi[512], al, bt, c[512] ;
  gdouble ui[3], vi[3], xi[453*4], ci[453*3], u, v, x[3], y[3] ;
  gdouble u1, u2, v1, v2, v3 ;
  gdouble emax ;
  
  fprintf(stderr, "spherical patch test\n") ;
  fprintf(stderr, "====================\n") ;
  fprintf(stderr, "nq = %d\n", nq) ;
  
  sqt_quadrature_select(nq, &q, &order) ;

  N = sqt_koornwinder_interp_matrix(q, nq, K) ;

  fprintf(stderr, "Knm N max: %d\n", N) ;

  u1 = 0.0 ; u2 = 0.7 ;
  v1 = 0.1 ; v2 = 0.1 ; v3 = 0.5 ;
  
  xstr = 4 ;
  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[i*3+0] ; t = q[i*3+1] ;

    u = u1*(1.0 - s - t) + u2*s + u2*t ;
    v = v1*(1.0 - s - t) + v2*s + v3*t ;

    sphere(u, v, &(xi[i*xstr])) ;
  }
  
#if 0

  /*compute the interpolation coefficients*/
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr, bt, ci, i3) ;

  emax = 0.0 ;
  for ( s = 0.0 ; s <= 1.0 ; s += 0.05 ) {
    for ( t = 0.0 ; t <= 1.0 - s ; t += 0.05 ) {
      u = ui[0]*(1.0 - s - t) + ui[1]*s + ui[2]*t ;
      v = vi[0]*(1.0 - s - t) + vi[1]*s + vi[2]*t ;
  
      stellarator(u, v, x, xu, xv) ;
      fprintf(stdout, "%lg %lg %lg ", x[0], x[1], x[2]) ;

      /*interpolate using K*/
      sqt_koornwinder_nm(N, s, t, Knm, 1) ;
      blaswrap_dgemv(TRUE, nq, i3, al, ci, i3, Knm, i1, bt, y, i1) ;
      
      fprintf(stdout, "%lg %lg %lg\n", y[0], y[1], y[2]) ;
      emax = MAX(emax,
		 (x[0]-y[0])*(x[0]-y[0]) +
		 (x[1]-y[1])*(x[1]-y[1]) +
		 (x[2]-y[2])*(x[2]-y[2])) ;
    }
  }

  fprintf(stderr, "maximum interpolation error: %lg\n", sqrt(emax)) ;
#endif
  
  return 0 ;
}

static gint blas_tests(gint N)

{
  gint i, j, stra, strx, stry, nr, nc, i1 = 1 ;
  gdouble A[8192], x[256], y[256], yref[256], al, bt, d, dref ;

  fprintf(stderr, "BLAS test\n") ;
  fprintf(stderr, "=========\n") ;

  nr = N ; nc = nr + 3 ;

  al = -1.23 ; bt = 0.7 ;
  
  stra = nc+5 ; strx = 2 ; stry = 5 ;
  /* stra = nc ; strx = 1 ; stry = 1 ; */
  
  fprintf(stderr, "A: %dx%d, stride %d\n", nr, nc, stra) ;
  fprintf(stderr, "x: %d elements, stride %d\n", nc, strx) ;
  fprintf(stderr, "y: %d elements, stride %d\n", nr, stry) ;
  fprintf(stderr, "al = %lg; bt = %lg\n", al, bt) ;
  
  for ( i = 0 ; i < nr ; i ++ ) {
    for ( j = 0 ; j < nc ; j ++ ) {
      A[i*stra+j] = (gdouble)(i+3)*(j+1) ;
    }
  }
  for ( j = 0 ; j < nc ; j ++ ) {
    y[stry*j] = -(gdouble)((j-0.3)*0.7) ;
    x[strx*j] = -(gdouble)((j+5.3)*0.7) ;
  }

  for ( i = 0 ; i < nr ; i ++ ) {
    yref[i] = bt*y[i*stry] ;
    for ( j = 0 ; j < nc ; j ++ ) {
      yref[i] += al*A[i*nc+j]*x[j*strx] ;
    }
  }

  blaswrap_dgemv(FALSE, nr, nc, al, A, stra, x, strx, bt, y, stry) ;

  for ( i = 0 ; i < nr ; i ++ ) {
    fprintf(stderr, "%lg ", yref[i]) ;
  }
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nr ; i ++ ) {
    fprintf(stderr, "%lg ", y[stry*i]) ;
  }
  fprintf(stderr, "\n") ;
  dref = 0.0 ;
  for ( i = 0 ; i < nr ; i ++ ) {
    fprintf(stderr, "%lg ", fabs(y[stry*i]-yref[i])) ;
    dref += y[stry*i]*yref[i] ;
  }
  fprintf(stderr, "\n") ;

  d = blaswrap_ddot(nr, y, stry, yref, i1) ;

  fprintf(stderr, "dot: %lg %lg (%lg)\n", d, dref, fabs(d-dref)) ;
  
  return 0 ;
}

static gint adaptive_quad_func(gdouble s, gdouble t, gdouble w,
			       gdouble *y, gdouble *n,
			       gdouble *quad, gint nq, gpointer data[])
{
  gdouble *x = data[0] ;
  gdouble R, dR, L[32] ;
  gint i ;
  
  g_assert(nq == 6 || nq == 12) ;
  g_assert(x != NULL) ;
  
  R = sqt_vector_distance(x, y) ;

  dR = sqt_vector_diff_scalar(x, y, n)/R/R ;  
  w *= 0.25*M_1_PI/R ;

  SQT_FUNCTION_NAME(sqt_element_shape_3d)(nq/2, s, t, L,
					  NULL, NULL, NULL, NULL, NULL) ;

  for ( i = 0 ; i < nq/2 ; i ++ ) {
    quad[     i] += w*L[i] ;
    quad[nq/2+i] += w*L[i]*dR ;
  }
  
  return 0 ;
}

static gint laplace_quad_func(gdouble s, gdouble t, gdouble w,
			      gdouble *y, gdouble *n,
			      gdouble *quad, gint nq, gpointer data[])

{
  gdouble *x = data[0] ;
  gdouble R, dR, L[32], src ;

  /* quad[0] += w ; */
  /* quad[1] += w ; */
  /* return 0 ; */
  
  g_assert(nq == 6 || nq == 12) ;
  /* g_assert(nq == 2) ; */
  g_assert(x != NULL) ;
  
  R = sqt_vector_distance(x, y) ;

  dR = sqt_vector_diff_scalar(x, y, n)/R/R ;  
  w *= 0.25*M_1_PI/R ;

  SQT_FUNCTION_NAME(sqt_element_shape_3d)(nq/2, s, t, L,
					  NULL, NULL, NULL, NULL, NULL) ;
  src = s*t - 1.0 ;
  quad[0] += w*src ;
  quad[1] += w*src*dR ;
  
  return 0 ;
}

static gint adaptive_quad_test(gdouble *xe, gint xstr, gint ne,
			       gint nq, gint N,
			       gdouble *x,
			       gint depth, gdouble tol,
			       gint nx)

{
  gdouble *q, f[512], g[512], t, *qref, qbas[512] ;
  gint oq, nc, i, nref ;
  sqt_quadrature_func_t func = (sqt_quadrature_func_t)adaptive_quad_func ;
  gpointer data[4] ;

  nc = 2*ne ;
  
  fprintf(stderr, "adaptive quadrature test\n") ;
  fprintf(stderr, "========================\n") ;
  fprintf(stderr, "x = %lg %lg %lg\n", x[0], x[1], x[2]) ;
  fprintf(stderr, "tol = %lg\n", tol) ;
  fprintf(stderr, "nq = %d\n", nq) ;

  nref = 453 ;
  sqt_quadrature_select(nref, &qref, &oq) ;
  sqt_quadrature_select(nq, &q, &oq) ;

  data[0] = x ;
  fprintf(stderr, "starting integration, t=%lg\n",
	  t = g_timer_elapsed(timer, NULL)) ;  
  sqt_adaptive_quad_tri(xe, xstr, ne, q, nq, func,
			f, nc, tol, depth, data) ;
  fprintf(stderr, "integration completed, t=%lg (%lg)\n",
	  g_timer_elapsed(timer, NULL), g_timer_elapsed(timer,NULL) - t) ;
  newman_tri_shape(x, &(xe[xstr*0]), &(xe[xstr*1]), &(xe[xstr*2]), NULL, 0,
		   &(g[0]), &(g[ne])) ;
  sqt_basic_quad_tri(xe, xstr, ne, qref, nref, func, qbas, nc, data) ;
  
  for ( i = 0 ; i < 2*ne ; i ++ ) g[i] *= -1.0/4.0/M_PI ;

  fprintf(stderr, "single layer\n") ;
  fprintf(stderr, "  adaptive:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", f[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  Newman: ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", g[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  basic:  ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", qbas[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error 1:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", fabs(f[i]-g[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error 2:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", fabs(f[i]-qbas[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "double layer\n") ;
  fprintf(stderr, "  adaptive:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", f[ne+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  Newman: ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", g[ne+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  basic:  ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", qbas[ne+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error 1:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg",
					fabs(f[ne+i]-g[ne+i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error 2:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg",
					fabs(f[ne+i]-qbas[ne+i])) ;
  fprintf(stderr, "\n") ;
  
  return 0 ;
}

static gint normal_quad_test(gdouble *xe, gint xstr, gint ne,
			     gint nq, gint N, gint depth, gdouble tol,
			     gdouble s0, gdouble t0, gdouble umax, gint nu)

{
  gdouble *q, J, n[3], g[512], x0[3], x[3], u ;
  gint oq, nc, i, j ;
  sqt_quadrature_func_t func = (sqt_quadrature_func_t)adaptive_quad_func ;
  gpointer data[4] ;

  nc = 2*ne ;
  
  sqt_element_point_3d(xe, xstr, ne, s0, t0, x0, n, &J) ;

  fprintf(stderr, "normal quadrature test\n") ;
  fprintf(stderr, "=====================\n") ;
  fprintf(stderr, "x0  = %lg %lg %lg\n", x0[0], x0[1], x0[2]) ;
  fprintf(stderr, "tol = %lg\n", tol) ;
  fprintf(stderr, "nq  = %d\n", nq) ;
  fprintf(stderr, "N   = %d\n", N) ;

  /*singular integral on surface to start with*/
  memset(g, 0, ne*sizeof(gdouble)) ;
  sqt_element_shape_3d(ne, s0, t0, &(g[ne]), NULL, NULL, NULL, NULL, NULL) ;
  for ( i = 0 ; i < ne ; i ++ ) g[ne+i] *= 0.5 ;  
  data[0] = x0 ;
  sqt_singular_quad_tri(xe, xstr, ne, s0, t0, N, func, g, nc, data) ;

  fprintf(stdout, "0") ;
  for ( i = 0 ; i < nc ; i ++ ) {
    fprintf(stdout, " %1.16e", g[i]) ;
  }
  fprintf(stdout, "\n") ;
  
  sqt_quadrature_select(nq, &q, &oq) ;
  data[0] = x ;

  for ( j = 1 ; j < nu ; j ++ ) {
    u = umax*j/(nu-1) ;
    x[0] = x0[0] + n[0]*u ;
    x[1] = x0[1] + n[1]*u ;
    x[2] = x0[2] + n[2]*u ;
    memset(g, 0, nc*sizeof(gdouble)) ;
    sqt_adaptive_quad_tri(xe, xstr, ne, q, nq, func, g, nc, tol, depth, data) ;
    fprintf(stdout, "%e", u) ;
    for ( i = 0 ; i < nc ; i ++ ) fprintf(stdout, " %1.16e", g[i]) ;
    fprintf(stdout, "\n") ;
  }

  return 0 ;
}

static gint quad_weight_test(gdouble *xe, gint xstr, gint ne,
			     gint nq, gint N, gint depth, gdouble tol,
			     gdouble s0, gdouble t0, gdouble umax, gint nu)

{
  gdouble *q, J, n[3], g[512], s, t ;
  gdouble x0[3], x[3], u, Kq[453*453], w[453*2], src[453] ;
  gdouble fs, fd, wb[453*2], err ;
  gint oq, nc, i, j, nK, nqk ;
  gpointer data[4] ;
  sqt_quadrature_func_t func = (sqt_quadrature_func_t)laplace_quad_func ;
  
  nc = 2*ne ;
  sqt_element_point_3d(xe, xstr, ne, s0, t0, x0, n, &J) ;

  nqk = 85 ;
  
  fprintf(stderr, "quadrature weight test\n") ;
  fprintf(stderr, "=====================\n") ;
  fprintf(stderr, "x0  = %lg %lg %lg\n", x0[0], x0[1], x0[2]) ;
  fprintf(stderr, "tol = %lg\n", tol) ;
  fprintf(stderr, "nq  = %d\n", nq) ;
  fprintf(stderr, "nqk = %d\n", nqk) ;


  sqt_quadrature_select(nqk, &q, &oq) ;
  nK = sqt_koornwinder_interp_matrix(q, nqk, Kq) ;
  for ( i = 0 ; i < nqk ; i ++ ) {
    s = q[3*i+0] ; t = q[3*i+1] ;
    src[i] = s*t - 1.0 ;
  }
  
  /*singular integral on surface to start with*/
  memset(g, 0, ne*sizeof(gdouble)) ;
  sqt_element_shape_3d(ne, s0, t0, &(g[ne]), NULL, NULL, NULL, NULL, NULL) ;

  g[1] = s0*t0 - 1.0 ;

  g[1] *= 0.5 ;
  
  data[0] = x0 ;
  sqt_laplace_weights_tri_singular(xe, xstr, ne, Kq, nqk, nK, N, s0, t0,
				   w) ;
  sqt_singular_quad_tri(xe, xstr, ne, s0, t0, N, func, g, nc, data) ;

  fs = 0.0 ; u = 0.0 ;
  fd = s0*t0 - 1.0 ;
  fd *= 0.5 ;
  for ( i = 0 ; i < nqk ; i ++ ) {
    fs += w[i]*src[i] ; fd += w[nqk+i]*src[i] ;
  }
  
  fprintf(stdout, "%e %e %e %e %e\n", u, fs, fd,
	  fabs(fs-g[0]), fabs(fd-g[1])) ;
  /* fprintf(stdout, "%e %e %e %e %e\n", u, fs, g[0], fd, g[1]) ; */
  
  data[0] = x ;
  
  sqt_quadrature_select(nq, &q, &oq) ;
  
  for ( j = 1 ; j < nu ; j ++ ) {
    u = umax*j/(nu-1) ;
    x[0] = x0[0] + n[0]*u ;
    x[1] = x0[1] + n[1]*u ;
    x[2] = x0[2] + n[2]*u ;

    memset(w, 0, 2*nqk*sizeof(gdouble)) ;
    sqt_laplace_weights_tri_adaptive(xe, xstr, ne, q, nq,
				     Kq, nqk, nK, tol, depth, x, w) ;
    g[0] = g[1] = 0.0 ;
    sqt_adaptive_quad_tri(xe, xstr, ne, q, nq, func, g, nc, tol, depth, data) ;

    fs = fd = 0.0 ;
    for ( i = 0 ; i < nqk ; i ++ ) {
      fs += w[i]*src[i] ; fd += w[nqk+i]*src[i] ;
    }

    fprintf(stdout, "%e %e %e %e %e\n", u, fs, fd,
	    fabs(fs-g[0]), fabs(fd-g[1])) ;
  }

  /*use last set of adaptively-integrated weights to check basic*/
  memset(wb, 0, 2*nqk*sizeof(gdouble)) ;
  nq = 175 ;
  
  sqt_quadrature_select(nq, &q, &oq) ;  
  sqt_laplace_weights_tri_basic(xe, xstr, ne, q, nq, Kq, nqk, nK, x, wb) ;

  err = 0.0 ;
  for ( i = 0 ; i < 2*nqk ; i ++ )
    err = MAX(err, fabs(wb[i] - w[i])) ;

  fprintf(stderr, "weights using %d point quadrature at u = %lg\n", nq, u) ;
  fprintf(stderr, "maximum error: %lg\n", err) ;
  
  return 0 ;
}

static gint matrix_adaptive_test(gdouble *xse, gint xsstr, gint nse,
				 gdouble *xte, gint xtstr, gint nte,
				 gint nq, gint depth, gdouble tol)
  

{
  gdouble *q, J, n[3], s, t, err ;
  gdouble x[3], Kq[453*453], src[453], al, bt, *st ;
  gdouble f[32], Ast[453*453*2], Astb[453*453*2], *qs, fts[453], ftd[453] ;
  gint oq, nc, i, nK, nqk, one = 1, nqt, nqs, lda ;
  gpointer data[4] ;
  sqt_quadrature_func_t func = (sqt_quadrature_func_t)laplace_quad_func ;

  nc = 2*nse ;

  nqk = 85 ;
  nqt = 25 ;
  nqs = 126 ;
  
  fprintf(stderr, "interaction matrix test\n") ;
  fprintf(stderr, "=======================\n") ;
  fprintf(stderr, "tol = %lg\n", tol) ;
  fprintf(stderr, "nqs = %d\n", nqs) ;
  fprintf(stderr, "nqk = %d\n", nqk) ;
  fprintf(stderr, "nqt = %d\n", nqt) ;

  sqt_quadrature_select(nqt, &st, &oq) ;

  sqt_quadrature_select(nqk, &q, &oq) ;
  nK = sqt_koornwinder_interp_matrix(q, nqk, Kq) ;
  for ( i = 0 ; i < nqk ; i ++ ) {
    s = q[3*i+0] ; t = q[3*i+1] ;
    src[i] = s*t - 1.0 ;
  }

  sqt_quadrature_select(nqs, &qs, &oq) ;

  /*interaction matrix*/
  sqt_laplace_source_target_tri_adaptive(xse, xsstr, nse, qs, nqs,
					 Kq, nqk, nK, tol, depth,
					 xte, xtstr, nte,
					 &(st[0]), 3, &(st[1]), 3, nqt,
					 Ast) ;

  /*calculate single and double layer potentials on target element*/
  al = 1.0 ; bt = 0.0 ; lda = 2*nqk ;
  blaswrap_dgemv(FALSE, nqt, nqk, al, &(Ast[0*nqk]), lda, src, one, bt,
		 fts, one) ;
  blaswrap_dgemv(FALSE, nqt, nqk, al, &(Ast[1*nqk]), lda, src, one, bt,
		 ftd, one) ;

  data[0] = x ;
  for ( i = 0 ; i < nqt ; i ++ ) {
    sqt_element_point_3d(xte, xtstr, nte, st[3*i+0], st[3*i+1], x, n, &J) ;
    sqt_adaptive_quad_tri(xse, xsstr, nse, qs, nqs, func, f, nc,
			  tol, depth, data) ;
    fprintf(stderr, "%lg %lg %lg %lg\n", fts[i], f[0], ftd[i], f[1]) ;
  }
  
  /*use last set of adaptively-integrated weights to check basic
   quadrature method for well-separated elements*/  
  sqt_quadrature_select(nq, &q, &oq) ;  
  sqt_laplace_source_target_tri_basic(xse, xsstr, nse, q, nq,
				      Kq, nqk, nK, xte, xtstr, nte,
				      &(st[0]), 3, &(st[1]), 3, nqt,
				      Astb) ;

  err = 0.0 ;
  for ( i = 0 ; i < 2*nqk*nqt ; i ++ )
    err = MAX(err, fabs(Ast[i] - Astb[i])) ;

  fprintf(stderr, "weights using %d point quadrature\n", nq) ;
  fprintf(stderr, "maximum error: %lg\n", err) ;
  
  return 0 ;
}

static gint matrix_self_test(gdouble *xe, gint xstr, gint ne, gint N)

{
  gdouble J, n[3], s, t, x[3] ;
  gdouble Kq[455*456*4], src[453], al, bt, *st ;
  gdouble f[32], fts[453], ftd[453], Ast[16384*4] ;
  gint oq, nc, i, nK, nqk, one = 1, lda ;
  gpointer data[4] ;
  sqt_quadrature_func_t func = (sqt_quadrature_func_t)laplace_quad_func ;

  nc = 2*ne ;

  nqk = 175 ;
  
  fprintf(stderr, "self-interaction matrix test\n") ;
  fprintf(stderr, "============================\n") ;
  fprintf(stderr, "N   = %d\n", N) ;
  fprintf(stderr, "nqk = %d\n", nqk) ;

  sqt_quadrature_select(nqk, &st, &oq) ;
  nK = sqt_koornwinder_interp_matrix(st, nqk, Kq) ;
  for ( i = 0 ; i < nqk ; i ++ ) {
    s = st[3*i+0] ; t = st[3*i+1] ;
    src[i] = s*t - 1.0 ;
    /* src[i] = 1.0 ; */
  }

  /*interaction matrix*/
  data[0] = x ;
  sqt_laplace_source_target_tri_self(xe, xstr, ne, Kq, nqk, nK, N,
				     &(st[0]), 3, &(st[1]), 3, nqk,
				     Ast) ;
  /*calculate single and double layer potentials on element*/
  al = 1.0 ; bt = 0.0 ; lda = 2*nqk ;
  blaswrap_dgemv(FALSE, nqk, nqk, al, &(Ast[0*nqk]), lda, src, one, bt,
  		 fts, one) ;
  blaswrap_dgemv(FALSE, nqk, nqk, al, &(Ast[1*nqk]), lda, src, one, bt,
  		 ftd, one) ;

  data[0] = x ;
  for ( i = 0 ; i < nqk ; i ++ ) {
    f[0] = f[1] = 0.0 ;
    sqt_element_point_3d(xe, xstr, ne, st[3*i+0], st[3*i+1], x, n, &J) ;
    sqt_singular_quad_tri(xe, xstr, ne, st[3*i+0], st[3*i+1], N, func,
    			  f, nc, data) ;
    fprintf(stderr, "%lg %lg %lg %lg (%e,%e)\n",
    	    fts[i], f[0], ftd[i], f[1],
    	    fabs(fts[i] - f[0]), fabs(ftd[i] - f[1])) ;
  }
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  gdouble xe[256], tol, rc, s0, t0, x0[3], umax, xt[256] ;
  gint ne, nte, nq, depth, N, nx, xstr, xtstr, test, i ;
  FILE *input ;
  gchar ch, *progname ;

  timer = g_timer_new() ;

  test = -1 ;
  depth = 0 ; tol = 1e-6 ; N = 8 ; nq = 54 ; rc = 0.05 ; nx = 33 ;
  umax = 0.1 ;
  s0 = t0 = G_MAXDOUBLE ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  input = stdin ;

  while ( (ch = getopt(argc, argv, "d:e:n:N:q:r:s:T:t:u:w:")) != EOF ) {
    switch (ch) {
    default: g_assert_not_reached() ; break ;
    case 'd': depth = atoi(optarg) ; break ;
    case 'e': tol   = atof(optarg) ; break ;
    case 'n': nx    = atoi(optarg) ; break ;
    case 'N': N     = atoi(optarg) ; break ;
    case 'q': nq    = atoi(optarg) ; break ;
    case 'r': rc    = atof(optarg) ; break ;
    case 's': s0    = atof(optarg) ; break ;
    case 'T': test  = parse_test(optarg) ; break ;      
    case 't': t0    = atof(optarg) ; break ;
    case 'u': umax  = atof(optarg) ; break ;
    }
  }

  if ( s0 != G_MAXDOUBLE && t0 == G_MAXDOUBLE ) t0 = s0 ;
  if ( t0 != G_MAXDOUBLE && s0 == G_MAXDOUBLE ) s0 = t0 ;

  if ( test == -1 ) {
    fprintf(stderr, "%s: unrecognized test case\n", progname) ;
    
    return 0 ;
  }
  
  if ( test == 1 ) {
    koornwinder_test(N, s0, t0) ;

    return 0 ;
  }

  if ( test == 2 ) {
    koornwinder_orthogonality_test(N) ;

    return 0 ;
  }

  if ( test == 3 ) {
    blas_tests(N) ;

    return 0 ;
  }
  
  if ( test == 4 ) {
    koornwinder_interpolation_test(N, nq) ;

    return 0 ;
  }

  if ( test == 10 ) {
    element_interpolation_test(N, nq) ;

    return 0 ;
  }

  if ( test == 11 ) {
    spherical_patch_test(N, nq) ;

    return 0 ;
  }

  if ( test == 12 ) {
    koornwinder_derivatives_test(N, s0, t0) ;

    return 0 ;
  }
  
  read_element(input, xe, &ne, &xstr) ;
  fscanf(input, "%lg %lg %lg", &(x0[0]), &(x0[1]), &(x0[2])) ;
	 
  fprintf(stderr, "%s: %d node element\n", progname, ne) ;

  if ( test == 0 ) {
    element_closest_point_test(xe, xstr, ne, s0, t0, x0, rc) ;
    
    return 0 ;
  }

  if ( test == 5 ) {
    adaptive_quad_test(xe, xstr, ne, nq, N, x0, depth, tol, nx) ;

    return 0 ;
  }
  
  if ( test == 6 ) {
    normal_quad_test(xe, xstr, ne, nq, N, depth, tol, s0, t0, umax, nx) ;

    return 0 ;
  }
  
  if ( test == 7 ) {
    quad_weight_test(xe, xstr, ne, nq, N, depth, tol, s0, t0, umax, nx) ;

    return 0 ;
  }

  if ( test == 8 ) {
    xtstr = xstr + 1 ; nte = ne ;
    for ( i = 0 ; i < ne ; i ++ ) {
      xt[xtstr*i+0] = xe[xstr*i+0] + 0.7 ; 
      xt[xtstr*i+1] = xe[xstr*i+1] + 1.7 ; 
      xt[xtstr*i+2] = xe[xstr*i+2] + 0.7 ; 
    }

    matrix_adaptive_test(xe, xstr, ne, xt, xtstr, nte, nq, depth, tol) ;

    return 0 ;
  }

  if ( test == 9 ) {

    matrix_self_test(xe, xstr, ne, N) ;

    return 0 ;
  }
  
  return 0 ;
}

