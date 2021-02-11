/*
 * Autogenerated file, do not  edit
 * Wed Feb 10 16:32:46 GMT 2021
 * -------------------------------
 */


/**
 * @ingroup location
 *
 * @brief Locate nearest point on an element
 *
 * Locate the nearest point on an element using the iterative
 * algorithm of Li X, Wu Z, Pan F et al. A geometric strategy
 * algorithm for orthogonal projection onto a parametric
 * surface. https://dx.doi.org/10.1007/s11390-019-1967-z
 *
 * @param xe element node coordinates;
 * @param xstr stride between \f$x\f$ coordinates of element nodes;
 * @param ne number of element nodes;
 * @param xf field point;
 * @param sn on exit reference element coordinate of nearest point;
 * @param tn on exit reference element coordinate of nearest point;
 * @param p1 physical location of nearest point;
 * @param tol tolerance for iteration;
 * @param nimax maximum number of iterations;
 * @param kn curvature at nearest point.
 *
 * @return number of iterations
 **/

gint sqt_element_nearest_point(gdouble *xe, gint xstr,
						  gint ne,
						  gdouble *xf, 
						  gdouble *sn, gdouble *tn, 
						  gdouble *p1, gdouble tol,
						  gint nimax, gdouble *kn) ;

/**
 * @ingroup base
 *
 * @brief Shape function evaluations for three-dimensional elements
 * 
 * On entry \a L may not be NULL, but if \a dLds is NULL, first and
 * second derivatives are not evaluated. If \a dLdss is NULL, first
 * derivatives are evaluated, but not second.
 * 
 * @param ne number of nodes on element (currently 3 or 6);
 * @param s coordinate on element;
 * @param t coordinate on element;
 * @param L on exit shape function evaluated at each node;
 * @param dLds on exit \f$\partial L_{i}/\partial s\f$ at each node;
 * @param dLdt on exit \f$\partial L_{i}/\partial t\f$ at each node;
 * @param dLdss on exit \f$\partial^{2} L_{i}/\partial s^{2}\f$ at each node;
 * @param dLdst on exit \f$\partial^{2} L_{i}/\partial s\partial t\f$ 
 * at each node;
 * @param dLdtt on exit \f$\partial^{2} L_{i}/\partial t^{2}\f$ at each node;
 *
 * @return 0 on success.
 * 
 **/

gint sqt_element_shape_3d(gint ne, gdouble s, gdouble t,
					     gdouble *L,
					     gdouble *dLds, gdouble *dLdt,
					     gdouble *dLdss, gdouble *dLdst,
					     gdouble *dLdtt) ;

/**
 * @ingroup base
 *
 * @brief Interpolation of point on element using shape functions
 * 
 * @param xe element node coordinates;
 * @param xstr stride between \f$x\f$ coordinates of element nodes;
 * @param ne number of element nodes;
 * @param L shape function from ::sqt_element_shape_3d(...);
 * @param dLds derivative of shape function;
 * @param dLdt derivative of shape function;
 * @param y on exit interpolated point on element;
 * @param n on exit normal at \a y;
 * @param J on exit Jacobian at \a y;
 * 
 * @return 0 on success.
 *
 **/

gint sqt_element_point_interp_3d(gdouble *xe,
						    gint xstr, gint ne,
						    gdouble *L,
						    gdouble *dLds,
						    gdouble *dLdt,
						    gdouble *y,
						    gdouble *n,
						    gdouble *J) ;

/**
 * @ingroup base
 *
 * @brief Convenience function for direct evaluation of point position
 * on element
 * 
 * @param xe element node coordinates;
 * @param xstr stride between \f$x\f$ coordinates of element nodes;
 * @param ne number of element nodes;
 * @param s coordinate on element;
 * @param t coordinate on element;
 * @param y on exit interpolated point on element;
 * @param n on exit normal at \a y;
 * @param J on exit Jacobian at \a y;
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_element_point_3d(gdouble *xe, gint xstr, gint ne,
					     gdouble s, gdouble t,
					     gdouble *y, gdouble *n,
					     gdouble *J) ;

/**
 * @ingroup base
 *
 * @brief Estimate surface area of element
 * 
 * @param xe element node coordinates;
 * @param xstr stride between \f$x\f$ coordinates of element nodes;
 * @param ne number of element nodes;
 * @param qrule nodes and weights of quadrature rule;
 * @param nq number of nodes in \a qrule.
 * 
 * @return estimate of area of element.
 * 
 **/

gdouble sqt_element_area(gdouble *xe, gint xstr, gint ne,
					     gdouble *qrule, gint nq) ;




/**
 * @ingroup base
 *
 * @brief Curvatures on a triangular element
 * 
 * @param xe element node coordinates;
 * @param xstr stride between \f$x\f$ coordinates of element nodes;
 * @param ne number of element nodes;
 * @param s coordinate on reference element;
 * @param t coordinate on reference element;
 * @param kg on exit Gaussian curvature on element in physical space at 
 * \f$(s,t)\f$;
 * @param km on exit mean curvature on element in physical space at
 * \f$(s,t)\f$.
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_triangle_curvature(gdouble *xe, gint xstr, gint ne,
					       gdouble s, gdouble t,
					       gdouble *kg, gdouble *km) ;

/**
 * @ingroup util
 *
 * @brief Convert Cartesian to spherical coordinates
 * 
 * @param x0 centre of coordinate system;
 * @param x point;
 * @param r radius of \a x in spherical system centred on \a x0;
 * @param th azimuth of \a x in spherical system centred on \a x0;
 * @param ph elevation of \a x in spherical system centred on \a x0.
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_cartesian_to_spherical(gdouble *x0,
						   gdouble *x,
						   gdouble *r,
						   gdouble *th,
						   gdouble *ph) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of Laplace potential from
 * triangular element using standard quadrature
 *
 * 
 * @param xe element nodes as three packed real numbers
 * @param xstr stride between nodes
 * @param ne number of nodes 
 * @param q nodes and weights of a triangular quadrature rule
 * @param nq number of nodes in \a q
 * @param Kq Koornwinder interpolation matrix from 
 * sqt_koornwinder_interp_matrix(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param x field point
 * @param w on exit contains weights for Laplace potential at \a x
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_laplace_weights_tri_basic(gdouble *xe,
						      gint xstr, gint ne,
						      gdouble *q, gint nq,
						      gdouble *Kq,
						      gint nqk,
						      gint nK,
						      gdouble *x,
						      gdouble *w) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of Laplace potential from
 * triangular element using adaptive quadrature
 *
 * 
 * @param xe element nodes as three packed real numbers
 * @param xstr stride between nodes
 * @param ne number of nodes 
 * @param q nodes and weights of a triangular quadrature rule
 * @param nq number of nodes in \a q
 * @param Kq Koornwinder interpolation matrix from 
 * sqt_koornwinder_interp_matrix(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param tol convergence tolerance for adaptive quadrature
 * @param dmax maximum recursion depth for adaptive quadrature
 * @param x field point
 * @param w on exit contains weights for Laplace potential at \a x
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_laplace_weights_tri_adaptive(gdouble *xe,
							 gint xstr, gint ne,
							 gdouble *q, gint nq,
							 gdouble *Kq,
							 gint nqk,
							 gint nK,
							 gdouble tol,
							 gint dmax,
							 gdouble *x,
							 gdouble *w) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of self-induced Laplace
 * potential on triangular element using singular quadrature
 *
 * 
 * @param xe element nodes as three packed real numbers
 * @param xstr stride between nodes
 * @param ne number of nodes 
 * @param Kq Koornwinder interpolation matrix from 
 * sqt_koornwinder_interp_matrix(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param N order of integration rules
 * @param s0 local coordinate of evaluation point on \a xe
 * @param t0 local coordinate of evaluation point on \a xe
 * @param w on exit contains weights for Laplace potential at \a x
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_laplace_weights_tri_singular(gdouble *xe,
							 gint xstr, gint ne,
							 gdouble *Kq,
							 gint nqk,
							 gint nK,
							 gint N,
							 gdouble s0,
							 gdouble t0,
							 gdouble *w) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of Laplace potential on
 * target triangular element from source triangular element using
 * standard quadrature
 *
 * 
 * @param xse source element nodes as three packed real numbers
 * @param xsstr stride between source element nodes
 * @param nse number of source element nodes 
 * @param q nodes and weights of a triangular quadrature rule
 * @param nq number of nodes in \a q
 * @param Kq Koornwinder interpolation matrix for \a xse from 
 * sqt_koornwinder_interp_matrix(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param xte target element nodes as three packed real numbers
 * @param xtstr stride between target element nodes
 * @param nte number of target element nodes 
 * @param s local coordinates of evaluation points on \a xte
 * @param sstr stride in \a s
 * @param t local coordinates of evaluation points on \a xte
 * @param tstr stride in \a t
 * @param nt number of evaluation points on \a xte
 * @param Ast on exit contains matrix for evaluation of single and double
 * layer potentials generated by source element on target
 * 
 * @return 0 on success.
 * 
 **/

gint
sqt_laplace_source_target_tri_basic(gdouble *xse,
						       gint xsstr, gint nse,
						       gdouble *q, gint nq,
						       gdouble *Kq,
						       gint nqk, gint nK,
						       gdouble *xte,
						       gint xtstr, gint nte,
						       gdouble *s,
						       gint sstr,
						       gdouble *t,
						       gint tstr,
						       gint nt,
						       gdouble *Ast) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of Laplace potential on
 * target triangular element from source triangular element using
 * adaptive quadrature
 * 
 * @param xse source element nodes as three packed real numbers
 * @param xsstr stride between source element nodes
 * @param nse number of source element nodes 
 * @param q nodes and weights of a triangular quadrature rule
 * @param nq number of nodes in \a q
 * @param Kq Koornwinder interpolation matrix for \a xse from 
 * sqt_koornwinder_interp_matrix(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param tol convergence tolerance for adaptive quadrature
 * @param dmax maximum recursion depth
 * @param xte target element nodes as three packed real numbers
 * @param xtstr stride between target element nodes
 * @param nte number of target element nodes 
 * @param s local coordinates of evaluation points on \a xte
 * @param sstr stride in \a s
 * @param t local coordinates of evaluation points on \a xte
 * @param tstr stride in \a t
 * @param nt number of evaluation points on \a xte
 * @param Ast on exit contains matrix for evaluation of single and double
 * layer potentials generated by source element on target
 * 
 * @return 0 on success.
 * 
 **/

gint
sqt_laplace_source_target_tri_adaptive(gdouble *xse,
							  gint xsstr, gint nse,
							  gdouble *q, gint nq,
							  gdouble *Kq,
							  gint nqk, gint nK,
							  gdouble tol,
							  gint dmax,
							  gdouble *xte,
							  gint xtstr, gint nte,
							  gdouble *s,
							  gint sstr,
							  gdouble *t,
							  gint tstr,
							  gint nt,
							  gdouble *Ast) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of self-induced Laplace
 * potential on triangular element using singular integral rules
 * 
 * @param xe element nodes as three packed real numbers
 * @param xstr stride between element nodes
 * @param ne number of element nodes 
 * @param Kq Koornwinder interpolation matrix for \a xse from 
 * sqt_koornwinder_interp_matrix(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param N order of integration method
 * @param s local coordinates of evaluation points on \a xe
 * @param sstr stride in \a s
 * @param t local coordinates of evaluation points on \a xe
 * @param tstr stride in \a t
 * @param nt number of evaluation points on \a xe
 * @param Ast on exit contains matrix for evaluation of single and double
 * layer potentials on element
 * 
 * @return 0 on success.
 * 
 **/

gint
sqt_laplace_source_target_tri_self(gdouble *xe,
						      gint xstr, gint ne,
						      gdouble *Kq,
						      gint nqk, gint nK,
						      gint N,
						      gdouble *s,
						      gint sstr,
						      gdouble *t,
						      gint tstr,
						      gint nt,
						      gdouble *Ast) ;

/**
 * @ingroup
 *
 * @brief
 * 
 * @param 
 * 
 * @return 0 on success.
 * 
 **/

/**
 * @ingroup location
 *
 * @brief Locate nearest point on an element
 *
 * Locate the nearest point on an element using the iterative
 * algorithm of Li X, Wu Z, Pan F et al. A geometric strategy
 * algorithm for orthogonal projection onto a parametric
 * surface. https://dx.doi.org/10.1007/s11390-019-1967-z
 *
 * @param xe element node coordinates;
 * @param xstr stride between \f$x\f$ coordinates of element nodes;
 * @param ne number of element nodes;
 * @param xf field point;
 * @param sn on exit reference element coordinate of nearest point;
 * @param tn on exit reference element coordinate of nearest point;
 * @param p1 physical location of nearest point;
 * @param tol tolerance for iteration;
 * @param nimax maximum number of iterations;
 * @param kn curvature at nearest point.
 *
 * @return number of iterations
 **/

gint sqt_element_nearest_point_f(gfloat *xe, gint xstr,
						  gint ne,
						  gfloat *xf, 
						  gfloat *sn, gfloat *tn, 
						  gfloat *p1, gfloat tol,
						  gint nimax, gfloat *kn) ;

/**
 * @ingroup base
 *
 * @brief Shape function evaluations for three-dimensional elements
 * 
 * On entry \a L may not be NULL, but if \a dLds is NULL, first and
 * second derivatives are not evaluated. If \a dLdss is NULL, first
 * derivatives are evaluated, but not second.
 * 
 * @param ne number of nodes on element (currently 3 or 6);
 * @param s coordinate on element;
 * @param t coordinate on element;
 * @param L on exit shape function evaluated at each node;
 * @param dLds on exit \f$\partial L_{i}/\partial s\f$ at each node;
 * @param dLdt on exit \f$\partial L_{i}/\partial t\f$ at each node;
 * @param dLdss on exit \f$\partial^{2} L_{i}/\partial s^{2}\f$ at each node;
 * @param dLdst on exit \f$\partial^{2} L_{i}/\partial s\partial t\f$ 
 * at each node;
 * @param dLdtt on exit \f$\partial^{2} L_{i}/\partial t^{2}\f$ at each node;
 *
 * @return 0 on success.
 * 
 **/

gint sqt_element_shape_3d_f(gint ne, gfloat s, gfloat t,
					     gfloat *L,
					     gfloat *dLds, gfloat *dLdt,
					     gfloat *dLdss, gfloat *dLdst,
					     gfloat *dLdtt) ;

/**
 * @ingroup base
 *
 * @brief Interpolation of point on element using shape functions
 * 
 * @param xe element node coordinates;
 * @param xstr stride between \f$x\f$ coordinates of element nodes;
 * @param ne number of element nodes;
 * @param L shape function from ::sqt_element_shape_3d_f(...);
 * @param dLds derivative of shape function;
 * @param dLdt derivative of shape function;
 * @param y on exit interpolated point on element;
 * @param n on exit normal at \a y;
 * @param J on exit Jacobian at \a y;
 * 
 * @return 0 on success.
 *
 **/

gint sqt_element_point_interp_3d_f(gfloat *xe,
						    gint xstr, gint ne,
						    gfloat *L,
						    gfloat *dLds,
						    gfloat *dLdt,
						    gfloat *y,
						    gfloat *n,
						    gfloat *J) ;

/**
 * @ingroup base
 *
 * @brief Convenience function for direct evaluation of point position
 * on element
 * 
 * @param xe element node coordinates;
 * @param xstr stride between \f$x\f$ coordinates of element nodes;
 * @param ne number of element nodes;
 * @param s coordinate on element;
 * @param t coordinate on element;
 * @param y on exit interpolated point on element;
 * @param n on exit normal at \a y;
 * @param J on exit Jacobian at \a y;
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_element_point_3d_f(gfloat *xe, gint xstr, gint ne,
					     gfloat s, gfloat t,
					     gfloat *y, gfloat *n,
					     gfloat *J) ;

/**
 * @ingroup base
 *
 * @brief Estimate surface area of element
 * 
 * @param xe element node coordinates;
 * @param xstr stride between \f$x\f$ coordinates of element nodes;
 * @param ne number of element nodes;
 * @param qrule nodes and weights of quadrature rule;
 * @param nq number of nodes in \a qrule.
 * 
 * @return estimate of area of element.
 * 
 **/

gfloat sqt_element_area_f(gfloat *xe, gint xstr, gint ne,
					     gfloat *qrule, gint nq) ;




/**
 * @ingroup base
 *
 * @brief Curvatures on a triangular element
 * 
 * @param xe element node coordinates;
 * @param xstr stride between \f$x\f$ coordinates of element nodes;
 * @param ne number of element nodes;
 * @param s coordinate on reference element;
 * @param t coordinate on reference element;
 * @param kg on exit Gaussian curvature on element in physical space at 
 * \f$(s,t)\f$;
 * @param km on exit mean curvature on element in physical space at
 * \f$(s,t)\f$.
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_triangle_curvature_f(gfloat *xe, gint xstr, gint ne,
					       gfloat s, gfloat t,
					       gfloat *kg, gfloat *km) ;

/**
 * @ingroup util
 *
 * @brief Convert Cartesian to spherical coordinates
 * 
 * @param x0 centre of coordinate system;
 * @param x point;
 * @param r radius of \a x in spherical system centred on \a x0;
 * @param th azimuth of \a x in spherical system centred on \a x0;
 * @param ph elevation of \a x in spherical system centred on \a x0.
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_cartesian_to_spherical_f(gfloat *x0,
						   gfloat *x,
						   gfloat *r,
						   gfloat *th,
						   gfloat *ph) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of Laplace potential from
 * triangular element using standard quadrature
 *
 * 
 * @param xe element nodes as three packed real numbers
 * @param xstr stride between nodes
 * @param ne number of nodes 
 * @param q nodes and weights of a triangular quadrature rule
 * @param nq number of nodes in \a q
 * @param Kq Koornwinder interpolation matrix from 
 * sqt_koornwinder_interp_matrix_f(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param x field point
 * @param w on exit contains weights for Laplace potential at \a x
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_laplace_weights_tri_basic_f(gfloat *xe,
						      gint xstr, gint ne,
						      gfloat *q, gint nq,
						      gfloat *Kq,
						      gint nqk,
						      gint nK,
						      gfloat *x,
						      gfloat *w) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of Laplace potential from
 * triangular element using adaptive quadrature
 *
 * 
 * @param xe element nodes as three packed real numbers
 * @param xstr stride between nodes
 * @param ne number of nodes 
 * @param q nodes and weights of a triangular quadrature rule
 * @param nq number of nodes in \a q
 * @param Kq Koornwinder interpolation matrix from 
 * sqt_koornwinder_interp_matrix_f(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param tol convergence tolerance for adaptive quadrature
 * @param dmax maximum recursion depth for adaptive quadrature
 * @param x field point
 * @param w on exit contains weights for Laplace potential at \a x
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_laplace_weights_tri_adaptive_f(gfloat *xe,
							 gint xstr, gint ne,
							 gfloat *q, gint nq,
							 gfloat *Kq,
							 gint nqk,
							 gint nK,
							 gfloat tol,
							 gint dmax,
							 gfloat *x,
							 gfloat *w) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of self-induced Laplace
 * potential on triangular element using singular quadrature
 *
 * 
 * @param xe element nodes as three packed real numbers
 * @param xstr stride between nodes
 * @param ne number of nodes 
 * @param Kq Koornwinder interpolation matrix from 
 * sqt_koornwinder_interp_matrix_f(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param N order of integration rules
 * @param s0 local coordinate of evaluation point on \a xe
 * @param t0 local coordinate of evaluation point on \a xe
 * @param w on exit contains weights for Laplace potential at \a x
 * 
 * @return 0 on success.
 * 
 **/

gint sqt_laplace_weights_tri_singular_f(gfloat *xe,
							 gint xstr, gint ne,
							 gfloat *Kq,
							 gint nqk,
							 gint nK,
							 gint N,
							 gfloat s0,
							 gfloat t0,
							 gfloat *w) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of Laplace potential on
 * target triangular element from source triangular element using
 * standard quadrature
 *
 * 
 * @param xse source element nodes as three packed real numbers
 * @param xsstr stride between source element nodes
 * @param nse number of source element nodes 
 * @param q nodes and weights of a triangular quadrature rule
 * @param nq number of nodes in \a q
 * @param Kq Koornwinder interpolation matrix for \a xse from 
 * sqt_koornwinder_interp_matrix_f(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param xte target element nodes as three packed real numbers
 * @param xtstr stride between target element nodes
 * @param nte number of target element nodes 
 * @param s local coordinates of evaluation points on \a xte
 * @param sstr stride in \a s
 * @param t local coordinates of evaluation points on \a xte
 * @param tstr stride in \a t
 * @param nt number of evaluation points on \a xte
 * @param Ast on exit contains matrix for evaluation of single and double
 * layer potentials generated by source element on target
 * 
 * @return 0 on success.
 * 
 **/

gint
sqt_laplace_source_target_tri_basic_f(gfloat *xse,
						       gint xsstr, gint nse,
						       gfloat *q, gint nq,
						       gfloat *Kq,
						       gint nqk, gint nK,
						       gfloat *xte,
						       gint xtstr, gint nte,
						       gfloat *s,
						       gint sstr,
						       gfloat *t,
						       gint tstr,
						       gint nt,
						       gfloat *Ast) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of Laplace potential on
 * target triangular element from source triangular element using
 * adaptive quadrature
 * 
 * @param xse source element nodes as three packed real numbers
 * @param xsstr stride between source element nodes
 * @param nse number of source element nodes 
 * @param q nodes and weights of a triangular quadrature rule
 * @param nq number of nodes in \a q
 * @param Kq Koornwinder interpolation matrix for \a xse from 
 * sqt_koornwinder_interp_matrix_f(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param tol convergence tolerance for adaptive quadrature
 * @param dmax maximum recursion depth
 * @param xte target element nodes as three packed real numbers
 * @param xtstr stride between target element nodes
 * @param nte number of target element nodes 
 * @param s local coordinates of evaluation points on \a xte
 * @param sstr stride in \a s
 * @param t local coordinates of evaluation points on \a xte
 * @param tstr stride in \a t
 * @param nt number of evaluation points on \a xte
 * @param Ast on exit contains matrix for evaluation of single and double
 * layer potentials generated by source element on target
 * 
 * @return 0 on success.
 * 
 **/

gint
sqt_laplace_source_target_tri_adaptive_f(gfloat *xse,
							  gint xsstr, gint nse,
							  gfloat *q, gint nq,
							  gfloat *Kq,
							  gint nqk, gint nK,
							  gfloat tol,
							  gint dmax,
							  gfloat *xte,
							  gint xtstr, gint nte,
							  gfloat *s,
							  gint sstr,
							  gfloat *t,
							  gint tstr,
							  gint nt,
							  gfloat *Ast) ;

/**
 * @ingroup laplace
 *
 * @brief Calculating weights for evaluation of self-induced Laplace
 * potential on triangular element using singular integral rules
 * 
 * @param xe element nodes as three packed real numbers
 * @param xstr stride between element nodes
 * @param ne number of element nodes 
 * @param Kq Koornwinder interpolation matrix for \a xse from 
 * sqt_koornwinder_interp_matrix_f(...)
 * @param nqk number of Koornwinder interpolation nodes
 * @param nK maximum order of Koornwinder polynomial
 * @param N order of integration method
 * @param s local coordinates of evaluation points on \a xe
 * @param sstr stride in \a s
 * @param t local coordinates of evaluation points on \a xe
 * @param tstr stride in \a t
 * @param nt number of evaluation points on \a xe
 * @param Ast on exit contains matrix for evaluation of single and double
 * layer potentials on element
 * 
 * @return 0 on success.
 * 
 **/

gint
sqt_laplace_source_target_tri_self_f(gfloat *xe,
						      gint xstr, gint ne,
						      gfloat *Kq,
						      gint nqk, gint nK,
						      gint N,
						      gfloat *s,
						      gint sstr,
						      gfloat *t,
						      gint tstr,
						      gint nt,
						      gfloat *Ast) ;

/**
 * @ingroup
 *
 * @brief
 * 
 * @param 
 * 
 * @return 0 on success.
 * 
 **/

