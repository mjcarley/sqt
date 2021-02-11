#ifndef SQT_H_INCLUDED
#define SQT_H_INCLUDED

#include <glib.h>

typedef gint (*sqt_quadrature_func_t)(gdouble s, gdouble t, gdouble w,
				    gdouble *y, gdouble *n,
				    gdouble *quad, gint nc, gpointer data) ;
typedef gint (*sqt_quadrature_func_f_t)(gfloat s, gfloat t, gfloat w,
				      gfloat *y, gfloat *n,
				      gfloat *quad, gint nc, gpointer data) ;
gint sqt_cartesian_to_spherical(gdouble *x0, gdouble *x, gdouble *r,
				gdouble *th, gdouble *ph) ;

gdouble sqt_minimum_distance(gdouble *xe, gint ne, gdouble *x) ;
gdouble sqt_maximum_distance(gdouble *xe, gint ne, gdouble *x) ;

gdouble sqt_element_area(gdouble *xe, gint xstr, gint ne,
			 gdouble *qrule, gint nq) ;
gint sqt_element_shape_3d(gint ne, gdouble s, gdouble t,
			  gdouble *L, gdouble *dLds, gdouble *dLdt,
			  gdouble *dLdss, gdouble *dLdst, gdouble *dLdtt) ;
gint sqt_element_point_3d(gdouble *xe, gint xstr, gint ne,
			  gdouble s, gdouble t,
			  gdouble *y, gdouble *n,
			  gdouble *J) ;
gint sqt_element_point_interp_3d(gdouble *xe, gint xstr, gint ne,
				 gdouble *L, gdouble *dLds, gdouble *dLdt,
				 gdouble *y, gdouble *n,
				 gdouble *J) ;

gint sqt_quadrature_select(gint nq, gdouble **q, gint *oq) ;

gint sqt_cartesian_to_spherical_f(gfloat *x0, gfloat *x, gfloat *r,
				  gfloat *th, gfloat *ph) ;
gfloat sqt_element_area_f(gfloat *xe, gint xstr, gint ne,
			  gfloat *qrule, gint nq) ;

gint sqt_element_shape_3d_f(gint ne, gfloat s, gfloat t,
			    gfloat *L, gfloat *dLds, gfloat *dLdt,
			    gfloat *dLdss, gfloat *dLdst, gfloat *dLdtt) ;
gint sqt_element_point_3d_f(gfloat *xe, gint xstr, gint ne,
			    gfloat s, gfloat t,
			    gfloat *y, gfloat *n,
			    gfloat *J) ;
gint sqt_element_point_interp_3d_f(gfloat *xe, gint xstr, gint ne,
				 gfloat *L, gfloat *dLds, gfloat *dLdt,
				 gfloat *y, gfloat *n,
				 gfloat *J) ;
gint sqt_triangle_curvature(gdouble *xe, gint xstr, gint ne,
			    gdouble s, gdouble t,
			    gdouble *k1, gdouble *k2) ;
gint sqt_triangle_curvature_f(gfloat *xe, gint xstr, gint ne,
			      gfloat s, gfloat t,
			      gfloat *k1, gfloat *k2) ;

gint sqt_quadrature_select_f(gint nq, gfloat **q, gint *oq) ;
gint sqt_element_nearest_point(gdouble *xe, gint xstr,
			       gint ne,
			       gdouble *xf, 
			       gdouble *sn, gdouble *tn, 
			       gdouble *xn, gdouble tol,
			       gint nimax, gdouble *kn) ;
gint sqt_element_nearest_point_f(gfloat *xe, gint xstr,
				 gint ne,
				 gfloat *xf, 
				 gfloat *sn, gfloat *tn, 
				 gfloat *xn, gfloat tol,
				 gint nimax, gfloat *kn) ;
gint sqt_koornwinder_nm(gint N, gdouble u, gdouble v, double *Knm,
			gint str, gint nk) ;
gint sqt_koornwinder_nm_f(gint N, gfloat u, gfloat v, gfloat *Knm,
			  gint str, gint nk) ;

gint sqt_koornwinder_deriv_nm(gint N, gdouble u, gdouble v,
			      gdouble *K , gint kstr,
			      gdouble *Ku, gint ustr,
			      gdouble *Kv, gint vstr, gint nk) ;
gint sqt_koornwinder_deriv_nm_f(gint N, gfloat u, gfloat v,
				gfloat *K , gint kstr,
				gfloat *Ku, gint ustr,
				gfloat *Kv, gint vstr, gint nk) ;

gint sqt_koornwinder_interp_matrix(gdouble *s, gint sstr,
				   gdouble *t, gint tstr,
				   gdouble *w, gint wstr,
				   gint nst, gdouble *A) ;
gint sqt_koornwinder_interp_matrix_f(gfloat *s, gint sstr,
				     gfloat *t, gint tstr,
				     gfloat *w, gint wstr,
				     gint nst, gfloat *A) ;

gint sqt_adaptive_quad_tri(gdouble *xe, gint xstr, gint ne,
			   gdouble *q, gint nq,
			   sqt_quadrature_func_t func,
			   gdouble *quad, gint nc, gdouble tol, gint dmax,
			   gpointer data) ;
gint sqt_adaptive_quad_tri_f(gfloat *xe, gint xstr, gint ne,
			     gfloat *q, gint nq,
			     sqt_quadrature_func_f_t func,
			     gfloat *quad, gint nc, gfloat tol, gint dmax,
			     gpointer data) ;

gint sqt_adaptive_quad_kw(gdouble *ce, gint ne, gint Nk,
			  gdouble *q, gint nq,
			  sqt_quadrature_func_t func,
			  gdouble *quad, gint nc,
			  gdouble tol, gint dmax,
			  gpointer data) ;
gint sqt_adaptive_quad_kw_f(gfloat *ce, gint ne, gint Nk,
			    gfloat *q, gint nq,
			    sqt_quadrature_func_f_t func,
			    gfloat *quad, gint nc,
			    gfloat tol, gint dmax,
			    gpointer data) ;
gint sqt_singular_quad_kw(gdouble *ce, gint ne, gint Nk,
			  gdouble s0, gdouble t0,
			  gint N,
			  sqt_quadrature_func_t func,
			  gdouble *quad, gint nc,
			  gpointer data) ;
gint sqt_singular_quad_kw_f(gfloat *ce, gint ne, gint Nk,
			    gfloat s0, gfloat t0,
			    gint N,
			    sqt_quadrature_func_f_t func,
			    gfloat *quad, gint nc,
			    gpointer data) ;

gint sqt_singular_quad_tri(gdouble *xe, gint xstr, gint ne,
			   gdouble s0, gdouble t0,
			   gint N,
			   sqt_quadrature_func_t func,
			   gdouble *quad, gint nc,
			   gpointer data) ;
gint sqt_singular_quad_tri_f(gfloat *xe, gint xstr, gint ne,
			     gfloat s0, gfloat t0,
			     gint N,
			     sqt_quadrature_func_f_t func,
			     gfloat *quad, gint nc,
			     gpointer data) ;

gint sqt_basic_quad_tri(gdouble *xe, gint xstr, gint ne,
			gdouble *q, gint nq,
			sqt_quadrature_func_t func,
			gdouble *quad, gint nc,
			gpointer data) ;

gint sqt_basic_quad_tri_f(gfloat *xe, gint xstr, gint ne,
			  gfloat *q, gint nq,
			  sqt_quadrature_func_f_t func,
			  gfloat *quad, gint nc,
			  gpointer data) ;

gint sqt_laplace_weights_tri_basic(gdouble *xe, gint xstr, gint ne,
				   gdouble *q, gint nq,
				   gdouble *Kq, gint nqk, gint nK,
				   gdouble *x, gdouble *w) ;
gint sqt_laplace_weights_tri_basic_f(gfloat *xe, gint xstr, gint ne,
				     gfloat *q, gint nq,
				     gfloat *Kq, gint nqk, gint nK,
				     gfloat *x, gfloat *w) ;

gint sqt_laplace_weights_tri_adaptive(gdouble *xe, gint xstr, gint ne,
				      gdouble *q, gint nq,
				      gdouble *Kq, gint nkq, gint nK,
				      gdouble tol, gint dmax,
				      gdouble *x,
				      gdouble *ws) ;
gint sqt_laplace_weights_tri_adaptive_f(gfloat *xe, gint xstr, gint ne,
					gfloat *q, gint nq,
					gfloat *Kq, gint nkq, gint nK,
					gfloat tol, gint dmax,
					gfloat *x,
					gfloat *ws) ;
gint sqt_laplace_weights_kw_adaptive(gdouble *ce, gint ne, gint Nk,
				     gdouble *K,
				     gdouble *q, gint nq,
				     gdouble tol, gint dmax,
				     gdouble *x,
				     gdouble *w) ;
gint sqt_laplace_weights_kw_adaptive_f(gfloat *ce, gint ne, gint Nk,
				       gfloat *K,
				       gfloat *q, gint nq,
				       gfloat tol, gint dmax,
				       gfloat *x,
				       gfloat *w) ;

gint sqt_laplace_weights_tri_singular(gdouble *xe, gint xstr, gint ne,
				      gdouble *Kq, gint nkq, gint nK,
				      gint N,
				      gdouble s0, gdouble t0,
				      gdouble *ws) ;
gint sqt_laplace_weights_tri_singular_f(gfloat *xe, gint xstr, gint ne,
					gfloat *Kq, gint nkq, gint nK,
					gint N,
					gfloat s0, gfloat t0,
					gfloat *gint) ;
gint sqt_laplace_weights_kw_singular(gdouble *ce, gint ne, gint Nk,
				     gdouble *Kq,
				     gint N,
				     gdouble s0,
				     gdouble t0,
				     gdouble *w) ;
gint sqt_laplace_weights_kw_singular_f(gfloat *ce, gint ne, gint Nk,
				       gfloat *Kq,
				       gint N,
				       gfloat s0,
				       gfloat t0,
				       gfloat *w) ;


gint sqt_laplace_source_target_tri_basic(gdouble *xse, gint xsstr, gint nse,
					 gdouble *q, gint nq,
					 gdouble *Kq, gint nqk, gint nK,
					 gdouble *xte, gint xtstr, gint nte,
					 gdouble *s, gint sstr,
					 gdouble *t, gint tstr, gint nt,
					 gdouble *Ast) ;
gint sqt_laplace_source_target_tri_basic_f(gfloat *xse, gint xsstr, gint nse,
					   gfloat *q, gint nq,
					   gfloat *Kq, gint nqk, gint nK,
					   gfloat *xte, gint xtstr, gint nte,
					   gfloat *s, gint sstr,
					   gfloat *t, gint tstr, gint nt,
					   gfloat *Ast) ;

gint sqt_laplace_source_target_tri_adaptive(gdouble *xse,
					    gint xsstr, gint nse,
					    gdouble *q, gint nq,
					    gdouble *Kq,
					    gint nqk, gint nK,
					    gdouble tol,
					    gint dmax,
					    gdouble *xte,
					    gint xtstr, gint nte,
					    gdouble *s, gint sstr,
					    gdouble *t, gint tstr,
					    gint nt,
					    gdouble *Ast) ;
gint sqt_laplace_source_target_tri_adaptive_f(gfloat *xse,
					      gint xsstr, gint nse,
					      gfloat *q, gint nq,
					      gfloat *Kq,
					      gint nqk, gint nK,
					      gfloat tol,
					      gint dmax,
					      gfloat *xte,
					      gint xtstr, gint nte,
					      gfloat *s, gint sstr,
					      gfloat *t, gint tstr,
					      gint nt,
					      gfloat *Ast) ;

gint sqt_laplace_source_target_kw_adaptive(gdouble *xse,
					   gint sstr, gint nse,
					   gdouble *q, gint nq,
					   gdouble *Ks, gint Ns,
					   gdouble tol,
					   gint dmax,
					   gdouble *xte,
					   gint tstr, gint nte,
					   gdouble *Ast) ;
gint sqt_laplace_source_target_kw_adaptive_f(gfloat *xse,
					     gint sstr, gint nse,
					     gfloat *q, gint nq,
					     gfloat *Ks, gint Ns,
					     gfloat tol,
					     gint dmax,
					     gfloat *xte,
					     gint tstr, gint nte,
					     gfloat *Ast) ;

gint sqt_laplace_source_target_tri_self(gdouble *xe, gint xstr, gint ne,
					gdouble *Kq, gint nqk, gint nK,
					gint N,
					gdouble *s, gint sstr,
					gdouble *t, gint tstr,
					gint nt,
					gdouble *Ast) ;
gint sqt_laplace_source_target_tri_self_f(gfloat *xe, gint xstr, gint ne,
					  gfloat *Kq, gint nqk, gint nK,
					  gint N,
					  gfloat *s, gint sstr,
					  gfloat *t, gint tstr,
					  gint nt,
					  gfloat *Ast) ;
gint sqt_laplace_source_target_kw_self(gdouble *xe, gint xstr, gint ne,
				       gdouble *K,
				       gint nK,
				       gint N,
				       gdouble *s,
				       gint sstr,
				       gdouble *t,
				       gint tstr,
				       gdouble *Ast) ;
gint sqt_laplace_source_target_kw_self_f(gfloat *xe, gint xstr, gint ne,
					 gfloat *K,
					 gint nK,
					 gint N,
					 gfloat *s,
					 gint sstr,
					 gfloat *t,
					 gint tstr,
					 gfloat *Ast) ;

gint sqt_element_interp(gdouble *ci, gint nq, gint Nk,
			gdouble s, gdouble t,
			gdouble *x, gdouble *n,
			gdouble *J, gdouble *dx, gdouble *work) ;
gint sqt_element_interp_f(gfloat *ci, gint nq, gint Nk,
			  gfloat s, gfloat t,
			  gfloat *x, gfloat *n,
			  gfloat *J, gfloat *dx, gfloat *work) ;

gint sqt_geometry_stellarator(gdouble u, gdouble v, gdouble *x, gdouble *n) ;
gint sqt_geometry_stellarator_f(gfloat u, gfloat v, gfloat *x, gfloat *n) ;

gint sqt_geometry_sphere(gdouble u, gdouble v, gdouble *x,
			 gdouble *xu, gdouble *xv) ;
gint sqt_geometry_sphere_f(gfloat u, gfloat v, gfloat *x,
			   gfloat *xu, gfloat *xv) ;

gint sqt_patch_nodes_tri(gdouble *xe, gint xstr, gint ne,
			 gdouble *s, gint sstr,
			 gdouble *t, gint tstr,
			 gint nst,
			 gdouble *xp, gint pstr,
			 gdouble *np, gint nstr) ;
gint sqt_patch_nodes_tri_f(gfloat *xe, gint xstr, gint ne,
			   gfloat *s, gint sstr,
			   gfloat *t, gint tstr,
			   gint nst,
			   gfloat *xp, gint pstr,
			   gfloat *np, gint nstr) ;
gint sqt_patch_nodes_sphere(gdouble rho,
			    gdouble th0, gdouble ph0,
			    gdouble th1, gdouble ph1,
			    gdouble th2, gdouble ph2,
			    gdouble *s, gint sstr,
			    gdouble *t, gint tstr,
			    gint nst,
			    gdouble *xp, gint pstr,
			    gdouble *np, gint nstr) ;
gint sqt_patch_nodes_sphere_f(gfloat rho,
			      gfloat th0, gfloat ph0,
			      gfloat th1, gfloat ph1,
			      gfloat th2, gfloat ph2,
			      gfloat *s, gint sstr,
			      gfloat *t, gint tstr,
			      gint nst,
			      gfloat *xp, gint pstr,
			      gfloat *np, gint nstr) ;

#endif /*SQT_H_INCLUDED*/
