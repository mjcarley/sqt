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
		  "matrix_indexed",
		  "cache",
		  "upsample",
		  "koornwinder_vector",
		  "koornwinder_deriv_vector",
		  "element_interpolation_vector",
		  ""} ;

GTimer *timer ;

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
  gint n, m, str, nq ;

  fprintf(stderr, "koornwinder test\n") ;
  fprintf(stderr, "================\n") ;

  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "(u,v) = (%lg,%lg)\n", u, v) ;

  str = 3 ; nq = 333 ;
  
  sqt_koornwinder_nm(N, u, v, Knm, str, nq) ;

  for ( n = 0 ; n < N ; n ++ ) {
    for ( m = 0 ; m <= n ; m ++ ) {
      if ( n*(n+1)/2+m < nq ) {
	fprintf(stdout,
		"%d %d %1.16e %1.16e %1.16e\n",
		n, m, u, v, Knm[str*(n*(n+1)/2+m)]) ;
      }
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
  str = 1 ;
  
  sqt_koornwinder_nm(N, u, v, K0, str, 1024) ;

  sqt_koornwinder_nm(N, u+0.5*ee, v       , Kup1, str, 1024) ;
  sqt_koornwinder_nm(N, u-0.5*ee, v       , Kum1, str, 1024) ;
  sqt_koornwinder_nm(N, u       , v+0.5*ee, Kvp1, str, 1024) ;
  sqt_koornwinder_nm(N, u       , v-0.5*ee, Kvm1, str, 1024) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    for ( m = 0 ; m <= n ; m ++ ) {
      idx = n*(n+1)/2 + m ;
      Kup1[str*idx] = (Kup1[str*idx] - Kum1[str*idx])/ee ;
      Kvp1[str*idx] = (Kvp1[str*idx] - Kvm1[str*idx])/ee ;
    }
  }
  
  sqt_koornwinder_deriv_nm(N, u, v, K, kstr, Ku, ustr, Kv, vstr, 1024) ;

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
	    
	    sqt_koornwinder_nm(N, s, t, Knm, str, nq) ;
	    
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
  gdouble s, t, Knm[453], K[453*453], *q, f, fr, fi[512], al, bt, c[512],
    emax ;
  
  fprintf(stderr, "koornwinder interpolation test\n") ;
  fprintf(stderr, "==============================\n") ;
  fprintf(stderr, "nq = %d\n", nq) ;
  
  sqt_quadrature_select(nq, &q, &order) ;

  N = sqt_koornwinder_interp_matrix(&(q[0]), 3, &(q[1]), 3, &(q[2]), 3, nq, K) ;

  fprintf(stderr, "Knm N max: %d\n", N) ;
  
  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[i*3+0] ; t = q[i*3+1] ;
    /* fi[i] = 3.0*s*t - t*t ; */
    fi[i] = sin(2.0*M_PI*s*t/8) ;
  }

  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemv(FALSE, nq, nq, al, K, nq, fi, i1, bt, c, i1) ;

  emax = 0.0 ;
  for ( s = 0 ; s <= 1.0 ; s += 0.1 ) {
    for ( t = 0 ; t <= 1.0-s ; t += 0.1 ) {
      sqt_koornwinder_nm(N, s, t, Knm, 1, nq) ;

      f = blaswrap_ddot(nq, c, i1, Knm, i1) ;
      /* fr = 3.0*s*t - t*t ; */
      fr = sin(2.0*M_PI*s*t/8) ;
  
      fprintf(stderr, "%lg %lg %lg %lg (%lg)\n", s, t, fr, f, fabs(fr-f)) ;
      emax = MAX(emax, fabs(fr-f)) ;
    }
  }

  fprintf(stderr, "maximum error: %lg\n", emax) ;
  
  return 0 ;
}

static gint element_interpolation_test(gint N, gint nq)

{
  gint order, i, i3 = 3, xstr ;
  gdouble s, t, K[4*65536], *q, al, bt ;
  gdouble ui[3], vi[3], xi[453*4], ci[453*3], u, v, x[3], xu[3], xv[3] ;
  gdouble y[3], n[3], ny[3], J ;
  gdouble emax, enmax ;
  gdouble work[3*453] ;
  
  fprintf(stderr, "element interpolation test\n") ;
  fprintf(stderr, "==========================\n") ;
  fprintf(stderr, "nq = %d\n", nq) ;
  
  sqt_quadrature_select(nq, &q, &order) ;

  N = sqt_koornwinder_interp_matrix(&(q[0]), 3, &(q[1]), 3, &(q[2]), 3, nq, K) ;

  fprintf(stderr, "Knm N max: %d\n", N) ;

  ui[0] = 0.3 ; ui[1] = 0.4  ; ui[2] = 0.35 ;
  vi[0] = 0.1 ; vi[1] = 0.11 ; vi[2] = 0.3 ;

  xstr = 4 ;
  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[i*3+0] ; t = q[i*3+1] ;

    u = ui[0]*(1.0 - s - t) + ui[1]*s + ui[2]*t ;
    v = vi[0]*(1.0 - s - t) + vi[1]*s + vi[2]*t ;

    sqt_geometry_stellarator(u, v, &(xi[i*xstr]), n) ;
  }

  /*compute the interpolation coefficients*/
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr, bt, ci, i3) ;

  emax = enmax = 0.0 ;
  for ( s = 0.01 ; s <= 1.0 ; s += 0.05 ) {
    for ( t = 0.01 ; t <= 1.0 - s ; t += 0.05 ) {
      u = ui[0]*(1.0 - s - t) + ui[1]*s + ui[2]*t ;
      v = vi[0]*(1.0 - s - t) + vi[1]*s + vi[2]*t ;
      
      sqt_geometry_stellarator(u, v, x, n) ;
      fprintf(stdout, "%lg %lg %lg\n", x[0], x[1], x[2]) ;

      sqt_element_interp(ci, nq, N, s, t, y, ny, &J, NULL, work) ;

      emax = MAX(emax,
		 (x[0]-y[0])*(x[0]-y[0]) +
		 (x[1]-y[1])*(x[1]-y[1]) +
		 (x[2]-y[2])*(x[2]-y[2])) ;
      enmax = MAX(enmax,
		  (n[0]-ny[0])*(n[0]-ny[0]) +
		  (n[1]-ny[1])*(n[1]-ny[1]) +
		  (n[2]-ny[2])*(n[2]-ny[2])) ;
    }
  }

  fprintf(stderr, "maximum location error: %lg\n", sqrt(emax)) ;
  fprintf(stderr, "maximum normal error: %lg\n", sqrt(enmax)) ;
  /* fprintf(stderr, "maximum differentiation error: %lg\n", sqrt(evmax)) ; */
  
  return 0 ;
}

static gint spherical_quad_func(gdouble s, gdouble t, gdouble w,
				gdouble *y, gdouble *n,
				gdouble *quad, gint nq, gint init,
				gpointer data[])
{
  gdouble R ;
  gdouble *x = data[0] ;
  
  R = sqt_vector_distance(x, y) ;

  quad[0] += w/R*0.25*M_1_PI ;
  quad[1] += w/R/R*0.25*M_1_PI ;
  
  return 0 ;
}

static gint spherical_test_func(gdouble th0, gdouble th1,
				gdouble ph0, gdouble ph1,
				gdouble rho,
				gint nc,
				gdouble *quad)

{
  gdouble a, gm, sc ;

  sc = 0.25*M_1_PI ;
  gm = (ph1*th1 - ph0*th0)/(th1 - th0) ;
  a = (ph1 - ph0)/(th1 - th0) ;

  quad[0] = (sin(gm - a*th1) - sin(gm - a*th0))/a
    + cos(ph0)*(th1 - th0) ;
  quad[0] *= sc*rho*rho/rho ;
  quad[1]  = quad[0]/rho ;
  
  return 0 ;
}

static gint spherical_patch_test(gint N, gint nq, gint depth, gdouble tol)

{
  gint order, i, i3 = 3, i1 = 1, xstr, nqk, nc, Nk ;
  gdouble K[454*175], *q, *qk, al, bt ;
  gdouble xi[175*4], ci[175*3], x[3], wt[175*2], src[175] ;
  gdouble th0, th1, ph0, ph1, rho, quad[1024], qref[1024], qwt[1024] ;
  gdouble *work ;
  sqt_quadrature_func_t func = (sqt_quadrature_func_t)spherical_quad_func ;
  gpointer data[8] ;
  
  fprintf(stderr, "spherical patch test\n") ;
  fprintf(stderr, "====================\n") ;
  fprintf(stderr, "nq = %d\n", nq) ;

  nqk = 175 ; nc = 2 ; xstr = 4 ;
  
  sqt_quadrature_select(nqk, &qk, &order) ;
  sqt_quadrature_select(nq, &q, &order) ;

  Nk = sqt_koornwinder_interp_matrix(&(qk[0]), 3, &(qk[1]), 3, &(qk[2]), 3,
				     nqk, K) ;

  fprintf(stderr, "Knm N max: %d\n", Nk) ;

  th0 = 0.4 ; th1 = 0.9 ;
  ph0 = 0.3 ; ph1 = 0.7 ;
  rho = 1.9 ;
  x[0] = x[1] = x[2] = 0.0 ;
  
  sqt_patch_nodes_sphere(rho, th0, ph0, th1, ph0, th0, ph1,
			 &(qk[0]), 3, &(qk[1]), 3, NULL, 1, nqk,
			 xi, xstr, NULL, 1, NULL, 1) ;

  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, nqk, i3, nqk, al, K, nqk, xi, xstr, bt, ci, i3) ;

  data[0] = x ;
  work = (gdouble *)g_malloc(4*2*nq*depth*sizeof(gdouble)) ;
  sqt_adaptive_quad_kw(ci, nqk, Nk, q, nq, func, quad, nc, tol, depth,
		       data, work) ;

  /*integration using pre-computed weights*/
  sqt_laplace_weights_kw_adaptive(ci, nqk, Nk, K, q, nq, tol, depth, x,
				  wt, work) ;
  
  for ( i = 0 ; i < nqk ; i ++ ) src[i] = 1.0 ;
  qwt[0] = blaswrap_ddot(nqk, &(wt[0  ]), i1, &(src[0]), i1) ;
  qwt[1] = blaswrap_ddot(nqk, &(wt[nqk]), i1, &(src[0]), i1) ;
  
  fprintf(stderr, "KW:      ") ;
  for ( i = 0 ; i < nc ; i ++ )
    fprintf(stderr, " %lg", quad[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "weights: ") ;
  for ( i = 0 ; i < nc ; i ++ )
    fprintf(stderr, " %lg", qwt[i]) ;
  fprintf(stderr, "\n") ;

  spherical_test_func(th0, th1, ph0, ph1, rho, nc, qref) ;
  fprintf(stderr, "ref:     ") ;
  for ( i = 0 ; i < nc ; i ++ )
    fprintf(stderr, " %lg", qref[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "err 1: ") ;
  for ( i = 0 ; i < nc ; i ++ )
    fprintf(stderr, " %lg", fabs(qref[i] - quad[i])) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "err 2: ") ;
  for ( i = 0 ; i < nc ; i ++ )
    fprintf(stderr, " %lg", fabs(qref[i] - qwt[i])) ;
  fprintf(stderr, "\n") ;

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
			       gdouble *quad, gint nq, gint init,
			       gpointer data[])
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
			      gdouble *quad, gint nq,
			      gint init, gpointer data[])

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

  /* SQT_FUNCTION_NAME(sqt_element_shape_3d)(nq/2, s, t, L, */
  /* 					  NULL, NULL, NULL, NULL, NULL) ; */
  src = s*t - 1.0 ;
  /* src = 1.0 ; */
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
  gdouble *q, f[512], g[512], *qref, qbas[512], qkw[512] ;
  gdouble xi[4096], ce[4096], K[453*453], al, bt, t ;
  gdouble *work ;
  gint oq, nc, i, nref, Nk, i3 = 3 ;
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

  /*Koornwinder interpolation data*/
  /* Nk = sqt_koornwinder_interp_matrix(q, nq, K) ; */
  Nk =
    sqt_koornwinder_interp_matrix(&(q[0]), 3, &(q[1]), 3, &(q[2]), 3, nq, K) ;
  sqt_patch_nodes_tri(xe, xstr, ne, &(q[0]), 3, &(q[1]), 3, NULL, 1, nq,
		      xi, xstr, NULL, 1, NULL, 1) ;
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr, bt, ce, i3) ;

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

  work = (gdouble *)g_malloc(4*nc*depth*sizeof(gdouble)) ;
  sqt_adaptive_quad_kw(ce, nq, Nk, q, nq, func, qkw, nc, tol, depth,
		       data, work) ;
  
  for ( i = 0 ; i < 2*ne ; i ++ ) g[i] *= -1.0/4.0/M_PI ;

  fprintf(stderr, "single layer\n") ;
  fprintf(stderr, "  adaptive:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", f[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  Newman:  ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", g[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  basic:   ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", qbas[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  KW   :   ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", qkw[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error 1:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", fabs(f[i]-g[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error 2:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", fabs(f[i]-qbas[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error 3:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", fabs(f[i]-qkw[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "double layer\n") ;
  fprintf(stderr, "  adaptive:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", f[ne+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  Newman:  ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", g[ne+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  basic:   ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", qbas[ne+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  KW   :   ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", qkw[ne+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error 1:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg",
					fabs(f[ne+i]-g[ne+i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error 2:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg",
					fabs(f[ne+i]-qbas[ne+i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error 2:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg",
					fabs(f[ne+i]-qkw[ne+i])) ;
  fprintf(stderr, "\n") ;
  
  return 0 ;
}

static gint cached_quad_test(gdouble *xe, gint xstr, gint ne,
			     gint nq, gint N,
			     gdouble *x,
			     gint depth, gdouble tol,
			     gint nx)

{
  gdouble *q, f[512], g[512], *qref, qbas[512], qkw[512], qkc[512] ;
  gdouble xi[4096], ce[4096], K[453*453], al, bt, t, *xcache ;
  gdouble *work ;
  gint oq, nc, i, nref, Nk, i3 = 3, *icache, cstr ;
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

  Nk =
    sqt_koornwinder_interp_matrix(&(q[0]), 3, &(q[1]), 3, &(q[2]), 3, nq, K) ;
  sqt_patch_nodes_tri(xe, xstr, ne, &(q[0]), 3, &(q[1]), 3, NULL, 1, nq,
		      xi, xstr, NULL, 1, NULL, 1) ;
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr, bt, ce, i3) ;

  data[0] = x ;

  cstr = NBI_CACHE_STRIDE + 7 ;
  icache = (gint *)g_malloc(nbi_cache_level_offset(depth+1)*sizeof(gint)) ;
  xcache = (gdouble *)g_malloc(nbi_cache_level_offset(depth+1)*nq*
			       cstr*sizeof(gdouble)) ;
  work = (gdouble *)g_malloc((4*nc*(depth)+3*nq)*sizeof(gdouble)) ;
  sqt_adaptive_quad_kw(ce, nq, Nk, q, nq, func, qkw, nc, tol, depth,
		       data, work) ;
  sqt_cached_quad_kw(ce, nq, Nk, q, nq, func, qkc, nc, tol, depth,
		     icache, xcache, cstr, TRUE, data, work) ;
  
  fprintf(stderr, "single layer\n") ;
  fprintf(stderr, "  cached:   ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", qkc[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  KW   :    ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", qkw[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", fabs(qkc[i]-qkw[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "double layer\n") ;

  fprintf(stderr, "  cached:   ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", qkc[ne+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  KW   :    ") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg", qkw[ne+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "  error:") ;
  for ( i = 0 ; i < ne ; i ++ ) fprintf(stderr, " %lg",
					fabs(qkc[ne+i]-qkw[ne+i])) ;
  fprintf(stderr, "\n") ;
  
  return 0 ;
}

static gint normal_quad_test(gdouble *xe, gint xstr, gint ne,
			     gint nq, gint N, gint depth, gdouble tol,
			     gdouble s0, gdouble t0, gdouble umax, gint nu)

{
  gdouble *q, J, n[3], g[512], gk[512], gw[512], K[453*453], ce[453*3] ;
  gdouble *qk, x0[3], x[3], u, xi[453*3], w[453*2] ;
  gdouble al, bt, ei, ew, eimax, ewmax ;
  gdouble src[6*453], *work ;
  gint oq, nc, i, j, i3 = 3, i1 = 1 ;
  gint nk, Nk, order ;
  sqt_quadrature_func_t func = (sqt_quadrature_func_t)adaptive_quad_func ;
  gpointer data[4] ;

  nc = 2*ne ;
  nk = 175 ;
    
  sqt_element_point_3d(xe, xstr, ne, s0, t0, x0, n, &J) ;

  fprintf(stderr, "normal quadrature test\n") ;
  fprintf(stderr, "=====================\n") ;
  fprintf(stderr, "x0   = %lg %lg %lg\n", x0[0], x0[1], x0[2]) ;
  fprintf(stderr, "tol  = %lg\n", tol) ;
  fprintf(stderr, "dmax = %d\n", depth) ;
  fprintf(stderr, "nq   = %d\n", nq) ;
  fprintf(stderr, "nk   = %d\n", nk) ;
  fprintf(stderr, "N    = %d\n", N) ;

  /*initialize Koornwinder interpolations*/
  sqt_quadrature_select(nk, &qk, &order) ;

  /* Nk = sqt_koornwinder_interp_matrix(qk, nk, K) ; */
  Nk = sqt_koornwinder_interp_matrix(&(qk[0]), 3, &(qk[1]), 3, &(qk[2]), 3,
				     nk, K) ;
  sqt_patch_nodes_tri(xe, xstr, ne, &(qk[0]), 3, &(qk[1]), 3, NULL, 1, nk,
		      xi, xstr, NULL, 1, NULL, 1) ;

  /*compute the interpolation coefficients*/
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, nk, i3, nk, al, K, nk, xi, xstr, bt, ce, i3) ;

  /*singular integral on surface to start with*/
  memset(g , 0, 2*ne*sizeof(gdouble)) ;
  memset(gk, 0, 2*ne*sizeof(gdouble)) ;
  memset(gw, 0, 2*ne*sizeof(gdouble)) ;
  sqt_element_shape_3d(ne, s0, t0, &(g[ne]), NULL, NULL, NULL, NULL, NULL) ;
  for ( i = 0 ; i < ne ; i ++ ) {
    g[ne+i] *= 0.5 ; gw[ne+i] = gk[ne+i] = g[ne+i] ;
  }
  data[0] = x0 ;
  sqt_singular_quad_tri(xe, xstr, ne, s0, t0, N, func, g, nc, data) ;
  sqt_singular_quad_kw(ce, nk, Nk, s0, t0, N, func, gk, nc, data) ;
  sqt_laplace_weights_kw_singular(ce, nk, Nk, K, N, s0, t0, w) ;
  for ( i = 0 ; i < nk ; i ++ ) {
    sqt_element_shape_3d(ne, qk[3*i+0], qk[3*i+1], &(src[i*ne]),
    			 NULL, NULL, NULL, NULL, NULL) ;
  }

  for ( i = 0 ; i < ne ; i ++ ) {
    gw[   i]  = blaswrap_ddot(nk, &(w[ 0]), i1, &(src[i]), ne) ;
    gw[ne+i] += blaswrap_ddot(nk, &(w[nk]), i1, &(src[i]), ne) ;
  }
  
  fprintf(stdout, "0") ; ei = ew = 0.0 ;
  for ( i = 0 ; i < nc ; i ++ ) {
    fprintf(stdout, " %lg", g[i]) ;
    ei = MAX(fabs(g[i] - gk[i]), ei) ;
    ew = MAX(fabs(g[i] - gw[i]), ew) ;
  }
  fprintf(stdout, " (%lg, %lg)\n", ei, ew) ;
  eimax = ei ; ewmax = ew ;
  
  sqt_quadrature_select(nq, &q, &oq) ;
  data[0] = x ;

  for ( j = 1 ; j < nu ; j ++ ) {
    u = umax*j/(nu-1) ;
    x[0] = x0[0] + n[0]*u ;
    x[1] = x0[1] + n[1]*u ;
    x[2] = x0[2] + n[2]*u ;
    memset(g , 0, nc*sizeof(gdouble)) ;
    memset(gk, 0, nc*sizeof(gdouble)) ;

    work = (gdouble *)g_malloc(4*2*ne*depth*sizeof(gdouble)) ;
    sqt_adaptive_quad_tri(xe, xstr, ne, q, nq, func, g, nc, tol, depth, data) ;
    sqt_adaptive_quad_kw(ce, nk, Nk, q, nq, func, gk, nc, tol, depth,
			 data, work) ;
    sqt_laplace_weights_kw_adaptive(ce, nk, Nk, K, q, nq, tol, depth, x,
				    w, work) ;
    for ( i = 0 ; i < ne ; i ++ ) {
      gw[   i] = blaswrap_ddot(nk, &(w[ 0]), i1, &(src[i]), ne) ;
      gw[ne+i] = blaswrap_ddot(nk, &(w[nk]), i1, &(src[i]), ne) ;
    }

    fprintf(stdout, "%lg", u) ; ei = ew = 0.0 ;
    for ( i = 0 ; i < nc ; i ++ ) {
      fprintf(stdout, " %lg", g[i]) ;
      ei = MAX(fabs(g[i] - gk[i]), ei) ;
      ew = MAX(fabs(g[i] - gw[i]), ew) ;
    }
    fprintf(stdout, " (%lg, %lg)\n", ei, ew) ;
    eimax = MAX(eimax, ei) ;
    ewmax = MAX(ewmax, ew) ;
  }

  fprintf(stderr, "maximum error: %lg, %lg\n", eimax, ewmax) ;
  
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
  nK = sqt_koornwinder_interp_matrix(&(q[0]), 3, &(q[1]), 3, &(q[2]), 3,
				     nqk, Kq) ;
  for ( i = 0 ; i < nqk ; i ++ ) {
    s = q[3*i+0] ; t = q[3*i+1] ;
    src[i] = s*t - 1.0 ;
  }
  
  /*singular integral on surface to start with*/
  memset(g, 0, ne*sizeof(gdouble)) ;
  sqt_element_shape_3d(ne, s0, t0, &(g[ne]), NULL, NULL, NULL, NULL, NULL) ;

  /* g[1] = s0*t0 - 1.0 ; */

  /* g[1] *= 0.5 ; */
  
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
  gdouble *q, J, n[3], s, t, err, time ;
  gdouble x[3], Kq[175*175], src[175], al, bt, *st ;
  gdouble Astb[175*175*2] ;
  gdouble f[32], Ast[175*175*2], xp[175*3], xt[175*3] ;
  gdouble *qs, fts[175], ftd[175], kts[175], ktd[175] ;
  gdouble *work ;
  gdouble ew, ek ;
  gint oq, nc, i, nK, nqk, one = 1, nqt, nqs, lda, pstr ;
  gpointer data[4] ;
  sqt_quadrature_func_t func = (sqt_quadrature_func_t)laplace_quad_func ;

  nc = 2*nse ;
  pstr = 3 ;
  
  nqk = 25 ;
  nqt = 54 ;
  nqs = 25 ;
  
  fprintf(stderr, "interaction matrix test\n") ;
  fprintf(stderr, "=======================\n") ;
  fprintf(stderr, "tol = %lg\n", tol) ;
  fprintf(stderr, "depth = %d\n", depth) ;
  fprintf(stderr, "nqs = %d\n", nqs) ;
  fprintf(stderr, "nqk = %d\n", nqk) ;
  fprintf(stderr, "nqt = %d\n", nqt) ;
  
  sqt_quadrature_select(nqt, &st, &oq) ;

  sqt_quadrature_select(nqk, &q, &oq) ;
  nK = sqt_koornwinder_interp_matrix(&(q[0]), 3, &(q[1]), 3, &(q[2]), 3,
				     nqk, Kq) ;
  for ( i = 0 ; i < nqk ; i ++ ) {
    s = q[3*i+0] ; t = q[3*i+1] ;
    src[i] = s*t - 1.0 ;
  }

  sqt_patch_nodes_tri(xse, xsstr, nse, &(q[0]), 3, &(q[1]), 3, NULL, 1, nqk,
  		      xp, pstr, NULL, 1, NULL, 1) ;
  sqt_patch_nodes_tri(xte, xtstr, nte, &(st[0]), 3, &(st[1]), 3, NULL, 1, nqt,
  		      xt, pstr, NULL, 1, NULL, 1) ;
  
  sqt_quadrature_select(nqs, &qs, &oq) ;

  /*interaction matrix*/
  sqt_laplace_source_target_tri_adaptive(xse, xsstr, nse, qs, nqs,
  					 Kq, nqk, nK, tol, depth,
  					 xte, xtstr, nte,
  					 &(st[0]), 3, &(st[1]), 3, nqt,
  					 Ast) ;
  work = (gdouble *)g_malloc(4*2*nqs*nqt*depth*sizeof(gdouble)) ;
  fprintf(stderr, "starting matrix generation, t=%lg\n",
	  time = g_timer_elapsed(timer, NULL)) ;    
  sqt_laplace_source_target_kw_adaptive(xp, pstr, nqk, qs, nqs, Kq, nK,
  					tol, depth, xt, pstr, nqt,
  					Astb, work) ;
  fprintf(stderr, "matrix generated, t=%lg\n",
	  g_timer_elapsed(timer, NULL) - time) ;    
  /* sqt_laplace_source_target_kw_adaptive(xp, pstr, nqk, qs, nqs, Kq, nK, */
  /* 					tol, depth, xt, pstr, nqt, */
  /* 					Astb) ; */
					
  /*calculate single and double layer potentials on target element*/
  al = 1.0 ; bt = 0.0 ; lda = 2*nqk ;
  blaswrap_dgemv(FALSE, nqt, nqk, al, &(Ast[0*nqk]), lda, src, one, bt,
		 fts, one) ;
  blaswrap_dgemv(FALSE, nqt, nqk, al, &(Ast[1*nqk]), lda, src, one, bt,
		 ftd, one) ;
  blaswrap_dgemv(FALSE, nqt, nqk, al, &(Astb[0*nqk]), lda, src, one, bt,
		 kts, one) ;
  blaswrap_dgemv(FALSE, nqt, nqk, al, &(Astb[1*nqk]), lda, src, one, bt,
		 ktd, one) ;

  data[0] = x ; ew = ek = 0.0 ;
  for ( i = 0 ; i < nqt ; i ++ ) {
    sqt_element_point_3d(xte, xtstr, nte, st[3*i+0], st[3*i+1], x, n, &J) ;
    sqt_adaptive_quad_tri(xse, xsstr, nse, qs, nqs, func, f, nc,
			  tol, depth, data) ;
    fprintf(stderr, "%lg %lg %lg %lg %lg %lg\n",
	    f[0], f[1],
	    fabs(fts[i]-f[0]), fabs(kts[i]-f[0]), 
	    fabs(ftd[i]-f[1]), fabs(ktd[i]-f[1])) ;
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
  fprintf(stderr, "maximum error in basic quadrature: %lg\n", err) ;

  return 0 ;
}

static gint matrix_indexed_test(gdouble *xse, gint xsstr, gint nse,
				gdouble *xte, gint xtstr, gint nte,
				gint nq, gint depth, gdouble tol)
  

{
  gdouble *q, err, erc, time ;
  gdouble Kq[453*453], *st ;
  gdouble *Ast, *Asti, *Astc, xp[453*3], xt[453*3], *xcache ;
  gdouble *qs, *work ;
  gint oq, i, j, k, nK, nqk, nqt, nqs, pstr, idx[453], ni, *icache, cstr ;

  pstr = 3 ;
  
  nqk = 7 ;
  nqt = 54 ;
  nqs = 25 ;

  Ast  = (gdouble *)g_malloc(2*nqk*nqt*sizeof(gdouble)) ;
  
  fprintf(stderr, "interaction matrix test\n") ;
  fprintf(stderr, "=======================\n") ;
  fprintf(stderr, "tol = %lg\n", tol) ;
  fprintf(stderr, "depth = %d\n", depth) ;
  fprintf(stderr, "nqs = %d\n", nqs) ;
  fprintf(stderr, "nqk = %d\n", nqk) ;
  fprintf(stderr, "nqt = %d\n", nqt) ;
  
  sqt_quadrature_select(nqt, &st, &oq) ;

  sqt_quadrature_select(nqk, &q, &oq) ;
  nK = sqt_koornwinder_interp_matrix(&(q[0]), 3, &(q[1]), 3, &(q[2]), 3,
				     nqk, Kq) ;
  sqt_patch_nodes_tri(xse, xsstr, nse, &(q[0]), 3, &(q[1]), 3, NULL, 1, nqk,
  		      xp, pstr, NULL, 1, NULL, 1) ;
  sqt_patch_nodes_tri(xte, xtstr, nte, &(st[0]), 3, &(st[1]), 3, NULL, 1, nqt,
  		      xt, pstr, NULL, 1, NULL, 1) ;
  
  sqt_quadrature_select(nqs, &qs, &oq) ;

  /*interaction matrix*/
  cstr = NBI_CACHE_STRIDE + nqs*2 ;
  work = (gdouble *)g_malloc(2*4*2*nqk*nqt*depth*sizeof(gdouble)) ;
  icache = (gint *)g_malloc(nbi_cache_level_offset(depth+1)*sizeof(gint)) ;
  xcache = (gdouble *)g_malloc(nbi_cache_level_offset(depth+1)*nqs*cstr*
			       sizeof(gdouble)) ;
  fprintf(stderr, "starting matrix generation, t=%lg\n",
	  time = g_timer_elapsed(timer, NULL)) ;    
  sqt_laplace_source_target_kw_adaptive(xp, pstr, nqk, qs, nqs, Kq, nK,
  					tol, depth, xt, pstr, nqt,
  					Ast, work) ;
  fprintf(stderr, "matrix generated, t=%lg\n",
	  g_timer_elapsed(timer, NULL) - time) ;    

  /*shuffle the indices*/
  ni = 2*nqt ;
  for ( i = 0 ; i < nqt ; i ++ ) {
    idx[i] = i ; idx[nqt+i] = nqt-i-1 ;
  }

  Asti = (gdouble *)g_malloc(2*nqk*ni*sizeof(gdouble)) ;
  Astc = (gdouble *)g_malloc(2*nqk*ni*sizeof(gdouble)) ;

  sqt_laplace_source_indexed_kw_adaptive(xp, pstr, nqk, qs, nqs, Kq, nK,
					 tol, depth, xt, pstr, idx, ni,
					 Asti, work) ;
  fprintf(stderr, "indexed matrix generated, t=%lg\n",
	  g_timer_elapsed(timer, NULL) - time) ;    
  sqt_laplace_source_indexed_kw_cached(xp, pstr, nqk, qs, nqs, Kq, nK,
				       tol, depth, xt, pstr, idx, ni,
				       icache, xcache, cstr,
				       Astc, work) ;
  fprintf(stderr, "indexed matrix generated (cached), t=%lg\n",
	  g_timer_elapsed(timer, NULL) - time) ;    

  err = erc = 0.0 ;
  for ( i = 0 ; i < ni ; i ++ ) {
    j = idx[i] ;
    for ( k = 0 ; k < 2*nqk ; k ++ ) {
      err = MAX(fabs(Ast[j*2*nqk+k] - Asti[i*2*nqk+k]), err) ;
      erc = MAX(fabs(Ast[j*2*nqk+k] - Astc[i*2*nqk+k]), erc) ;
    }
  }    
  
  fprintf(stderr, "%d indexed rows, error = %lg\n", ni, err) ;
  fprintf(stderr, "%d indexed rows, cached error = %lg\n", ni, erc) ;

  return 0 ;
}

static gint matrix_self_test(gdouble *xe, gint xstr, gint ne, gint N,
			     gint nqk)

{
  gdouble J, n[3], s, t, x[3] ;
  gdouble *Kq, *src, al, bt, *st, *xp ;
  gdouble f[32], *Ast, *Astk ;
  gdouble *kts, *ktd, *fts, *ftd ;
  gint oq, nc, i, nK, one = 1, lda, pstr ;
  gpointer data[4] ;
  sqt_quadrature_func_t func = (sqt_quadrature_func_t)laplace_quad_func ;

  nc = 2*ne ; pstr = 3 ;
  
  fprintf(stderr, "self-interaction matrix test\n") ;
  fprintf(stderr, "============================\n") ;
  fprintf(stderr, "N   = %d\n", N) ;
  fprintf(stderr, "nqk = %d\n", nqk) ;

  Kq   = (gdouble *)g_malloc(  nqk*nqk*sizeof(gdouble)) ;
  Ast  = (gdouble *)g_malloc(2*nqk*nqk*sizeof(gdouble)) ;
  Astk = (gdouble *)g_malloc(2*nqk*nqk*sizeof(gdouble)) ;

  fts  = (gdouble *)g_malloc(      nqk*sizeof(gdouble)) ;
  ftd  = (gdouble *)g_malloc(      nqk*sizeof(gdouble)) ;
  kts  = (gdouble *)g_malloc(      nqk*sizeof(gdouble)) ;
  ktd  = (gdouble *)g_malloc(      nqk*sizeof(gdouble)) ;
  src  = (gdouble *)g_malloc(      nqk*sizeof(gdouble)) ;
  xp   = (gdouble *)g_malloc(    3*nqk*sizeof(gdouble)) ;
  
  sqt_quadrature_select(nqk, &st, &oq) ;
  nK = sqt_koornwinder_interp_matrix(&(st[0]), 3, &(st[1]), 3, &(st[2]), 3,
				     nqk, Kq) ;
  for ( i = 0 ; i < nqk ; i ++ ) {
    s = st[3*i+0] ; t = st[3*i+1] ;
    src[i] = s*t - 1.0 ;
    /* src[i] = 1.0 ; */
  }
  sqt_patch_nodes_tri(xe, xstr, ne, &(st[0]), 3, &(st[1]), 3, NULL, 1, nqk,
  		      xp, pstr, NULL, 1, NULL, 1) ;
  sqt_laplace_source_target_kw_self(xp, pstr, nqk, Kq, nK, N,
  				    &(st[0]), 3, &(st[1]), 3, Astk) ;
  
  /*interaction matrix*/
  sqt_laplace_source_target_tri_self(xe, xstr, ne, Kq, nqk, nK, N,
				     &(st[0]), 3, &(st[1]), 3, nqk,
				     Ast) ;
  /*calculate single and double layer potentials on element*/
  al = 1.0 ; bt = 0.0 ; lda = 2*nqk ;
  blaswrap_dgemv(FALSE, nqk, nqk, al, &(Ast[0*nqk]), lda, src, one, bt,
  		 fts, one) ;
  blaswrap_dgemv(FALSE, nqk, nqk, al, &(Ast[1*nqk]), lda, src, one, bt,
  		 ftd, one) ;
  blaswrap_dgemv(FALSE, nqk, nqk, al, &(Astk[0*nqk]), lda, src, one, bt,
  		 kts, one) ;
  blaswrap_dgemv(FALSE, nqk, nqk, al, &(Astk[1*nqk]), lda, src, one, bt,
  		 ktd, one) ;

  data[0] = x ;
  for ( i = 0 ; i < nqk ; i ++ ) {
    f[0] = f[1] = 0.0 ;
    sqt_element_point_3d(xe, xstr, ne, st[3*i+0], st[3*i+1], x, n, &J) ;
    sqt_singular_quad_tri(xe, xstr, ne, st[3*i+0], st[3*i+1], N, func,
    			  f, nc, data) ;
    fprintf(stderr, "%lg %lg (%lg,%lg,%lg,%lg)\n",
    	    f[0], f[1],
    	    fabs(fts[i] - f[0]), fabs(ftd[i] - f[1]),
    	    fabs(kts[i] - f[0]), fabs(ktd[i] - f[1])) ;
    /* fprintf(stdout, "%lg %lg %lg %lg %lg %lg\n", */
    /* 	    f[0], f[1], kts[i], ktd[i], fts[i], ftd[i]) ; */
  }

  return 0 ;
}

gdouble upsample_test_func(gdouble s, gdouble t)

{
  gdouble f ;

  f = (1.0 - s*s)*t - 3.0*t*t*t ;
  
  return f ;
}
  
gint upsample_test(gint nq)

{
  gdouble *K, *st, *q, *Ki, f[453], g[453], s, t, work[453], al, bt, err ;
  gint order, Nk, nqi[] = {7, 25, 54, 85, 126, 175, 453}, i, j, i1 = 1 ;
  
  fprintf(stderr, "upsample matrix test\n") ;
  fprintf(stderr, "====================\n") ;
  fprintf(stderr, "nq = %d\n", nq) ;

  K  = (gdouble *)g_malloc0(nq*nq*sizeof(gdouble)) ;
  Ki = (gdouble *)g_malloc0(nq*453*sizeof(gdouble)) ;
  sqt_quadrature_select(nq, &q, &order) ;

  Nk = sqt_koornwinder_interp_matrix(&(q[0]), 3, &(q[1]), 3, &(q[2]), 3, nq,
				     K) ;

  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[3*i+0] ; t = q[3*i+1] ;
    f[i] = upsample_test_func(s, t) ;
  }

  al = 1.0 ; bt = 0.0 ;
  for ( i = 0 ; i < 7 ; i ++ ) {
    sqt_quadrature_select(nqi[i], &st, &order) ;
    sqt_interp_matrix(K, nq, Nk, &(st[0]), 3, &(st[1]), 3, nqi[i], Ki, work) ;
    blaswrap_dgemv(FALSE, nqi[i], nq, al, Ki, nq, f, i1, bt, g, i1) ;
    err = 0.0 ;
    for ( j = 0 ; j < nqi[i] ; j ++ ) {
      s = st[3*j+0] ; t = st[3*j+1] ; 
      err = MAX(err, fabs(g[j] - upsample_test_func(s, t))) ;
    }
    fprintf(stderr, "%d upsampled to %d nodes, error: %lg\n",
	    nq, nqi[i], err) ;
  }
  
  return 0 ;
}

static gint koornwinder_vector_test(gint N)

{
  gdouble Knm[1024], Kvec[4096], s[4], t[4], err ;
  gint n, m, str, nq, i, j, idx, nk ;

  fprintf(stderr, "koornwinder vector test\n") ;
  fprintf(stderr, "=======================\n") ;

  fprintf(stderr, "N = %d\n", N) ;

  str = 3 ; nq = 333 ; nk = 4 ;

  s[0] = 0.1 ; s[1] = 0.3 ; s[2] = 0.37 ; s[3] = 0.7 ;
  t[0] = 0.6 ; t[1] = 0.2 ; t[2] = 0.57 ; t[3] = 0.1 ;

  sqt_koornwinder_nm_vector(N, s, t, nk, Kvec, str, nq) ;

  for ( i = 0 ; i < nk ; i ++ ) {
    sqt_koornwinder_nm(N, s[i], t[i], Knm, str, nq) ;

    err = 0.0 ;
    for ( n = 0 ; n < N ; n ++ ) {
      for ( m = 0 ; m <= n ; m ++ ) {
	idx = n*(n+1)/2+m ;
	if ( idx < nq ) { 
	  err = MAX(err, fabs(Knm[idx*str] - Kvec[(i*nq+idx)*str])) ;
	}
      }
    }

    fprintf(stderr, "%d (%lg, %lg) %lg\n", i, s[i], t[i], err) ;
  }
  
  return 0 ;
}

static gint koornwinder_deriv_vector_test(gint N)

{
  gdouble Knm[1024], Ks[1024], Kt[1024], Kv[16384], Ksv[16384], Ktv[16384] ;
  gdouble s[4], t[4], errk, errs, errt ;
  gint n, m, str, nq, i, j, idx, nk, kstr, sstr, tstr, offk, offs, offt ;

  fprintf(stderr, "koornwinder deriv vector test\n") ;
  fprintf(stderr, "=============================\n") ;

  fprintf(stderr, "N = %d\n", N) ;

  nq = 333 ; nk = 4 ;
  kstr = 3 ; sstr = 2 ; tstr = 2 ;
  offk = nq+7 ; offs = nq+3 ; offt = nq+5 ;
  
  s[0] = 0.1 ; s[1] = 0.3 ; s[2] = 0.37 ; s[3] = 0.7 ;
  t[0] = 0.6 ; t[1] = 0.2 ; t[2] = 0.57 ; t[3] = 0.1 ;

  sqt_koornwinder_deriv_nm_vector(N, s, t, nk,
				  Kv, kstr, offk,
				  Ksv, sstr, offs,
				  Ktv, tstr, offt,
				  nq) ;

  for ( i = 0 ; i < nk ; i ++ ) {
    sqt_koornwinder_deriv_nm(N, s[i], t[i], Knm, kstr, Ks, sstr, Kt, tstr,
			     nq) ;
    /* sqt_koornwinder_nm(N, s[i], t[i], Knm, str, nq) ; */

    errs = errt = errk = 0.0 ;
    for ( n = 0 ; n < N ; n ++ ) {
      for ( m = 0 ; m <= n ; m ++ ) {
	idx = n*(n+1)/2+m ;
	if ( idx < nq ) { 
	  errk = MAX(errk, fabs(Knm[idx*kstr] - Kv[(i*offk+idx)*kstr])) ;
	  errs = MAX(errs, fabs(Ks[idx*sstr]  - Ksv[(i*offs+idx)*sstr])) ;
	  errt = MAX(errt, fabs(Kt[idx*tstr]  - Ktv[(i*offt+idx)*tstr])) ;
	}
      }
    }

    fprintf(stderr, "%d (%lg, %lg) %lg, %lg, %lg\n",
	    i, s[i], t[i], errk, errs, errt) ;
  }
  
  return 0 ;
}

static gint element_interpolation_vector_test(gint N, gint nq)

{
  gint order, i, i3 = 3, xstr ;
  gdouble s[4], t[4], K[4*65536], *q, al, bt, s0, t0 ;
  gdouble ui[3], vi[3], xi[453*4], ci[453*3], u, v, x[12], xu[12], xv[12] ;
  gdouble y[12], n[12], ny[12], J[4], Jy ;
  gdouble emax, enmax ;
  gdouble work[4*3*453] ;
  
  fprintf(stderr, "element interpolation vector test\n") ;
  fprintf(stderr, "=================================\n") ;
  fprintf(stderr, "nq = %d\n", nq) ;
  
  sqt_quadrature_select(nq, &q, &order) ;

  N = sqt_koornwinder_interp_matrix(&(q[0]), 3, &(q[1]), 3, &(q[2]), 3, nq, K) ;

  fprintf(stderr, "Knm N max: %d\n", N) ;

  ui[0] = 0.3 ; ui[1] = 0.4  ; ui[2] = 0.35 ;
  vi[0] = 0.1 ; vi[1] = 0.11 ; vi[2] = 0.3 ;

  xstr = 4 ;
  for ( i = 0 ; i < nq ; i ++ ) {
    s0 = q[i*3+0] ; t0 = q[i*3+1] ;

    u = ui[0]*(1.0 - s0 - t0) + ui[1]*s0 + ui[2]*t0 ;
    v = vi[0]*(1.0 - s0 - t0) + vi[1]*s0 + vi[2]*t0 ;

    sqt_geometry_stellarator(u, v, &(xi[i*xstr]), n) ;
  }

  /*compute the interpolation coefficients*/
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr, bt, ci, i3) ;

  s[0] = 0.1 ; s[1] = 0.3 ; s[2] = 0.37 ; s[3] = 0.7 ;
  t[0] = 0.6 ; t[1] = 0.2 ; t[2] = 0.57 ; t[3] = 0.1 ;

  sqt_element_interp_vector(ci, nq, N, s, t, 4, x, 3, n, 3, J, work) ;

  for ( i = 0 ; i < 4 ; i ++ ) {
    sqt_element_interp(ci, nq, N, s[i], t[i], y, ny, &Jy, NULL, work) ;
    fprintf(stderr, "%d: %lg, %lg, %lg, %lg, %lg, %lg (%lg)\n",
	    i, y[0], y[1], y[2], ny[0], ny[1], ny[2], Jy) ;
    fprintf(stderr, "   %lg, %lg, %lg, %lg, %lg, %lg (%lg)\n",
	    x[3*i+0], x[3*i+1], x[3*i+2], n[3*i+0], n[3*i+1], n[3*i+2], J[i]) ;
  }
    
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  gdouble xe[256], tol, rc, s0, t0, x0[3], umax, xt[256], xm[256] ;
  gint ne, nte, nq, depth, N, nx, xstr, xtstr, test, i ;
  FILE *input ;
  gchar ch, *progname ;

  timer = g_timer_new() ;

  memset(xe, 0, 256*sizeof(gdouble)) ;
  memset(xt, 0, 256*sizeof(gdouble)) ;
  
  test = -1 ;
  depth = 8 ; tol = 1e-6 ; N = 8 ; nq = 54 ; rc = 0.05 ; nx = 33 ;
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
    spherical_patch_test(N, nq, depth, tol) ;

    return 0 ;
  }

  if ( test == 12 ) {
    koornwinder_derivatives_test(N, s0, t0) ;

    return 0 ;
  }

  if ( test == 15 ) {
    upsample_test(nq) ;

    return 0 ;
  }

  if ( test == 16 ) {
    koornwinder_vector_test(N) ;

    return 0 ;
  }
  
  if ( test == 17 ) {
    koornwinder_deriv_vector_test(N) ;

    return 0 ;
  }
  
  if ( test == 18 ) {
    element_interpolation_vector_test(N, nq) ;

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

    matrix_self_test(xe, xstr, ne, N, nq) ;

    return 0 ;
  }
  
  if ( test == 13 ) {
    xtstr = xstr + 1 ; nte = ne ;
    for ( i = 0 ; i < ne ; i ++ ) {
      xt[xtstr*i+0] = xe[xstr*i+0] + 0.7 ; 
      xt[xtstr*i+1] = xe[xstr*i+1] + 1.7 ; 
      xt[xtstr*i+2] = xe[xstr*i+2] + 0.7 ; 
    }

    matrix_indexed_test(xe, xstr, ne, xt, xtstr, nte, nq, depth, tol) ;

    return 0 ;
  }

  if ( test == 14 ) {
    cached_quad_test(xe, xstr, ne, nq, N, x0, depth, tol, nx) ;
    
    return 0 ;
  }
  
  return 0 ;
}

