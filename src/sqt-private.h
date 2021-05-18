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

#ifndef SQT_PRIVATE_H_INCLUDED 
#define SQT_PRIVATE_H_INCLUDED

#include <stdio.h>

#include <glib.h>

#define SQT_DATA_WIDTH         16
#define SQT_DATA_ELEMENT        0
#define SQT_DATA_STRIDE         1
#define SQT_DATA_NUMBER         2
#define SQT_DATA_TARGET         3
#define SQT_DATA_NORMAL         4
#define SQT_DATA_MATRIX         5
#define SQT_DATA_KNM            6
#define SQT_DATA_NKNM           7
#define SQT_DATA_ORDER_K        8
#define SQT_DATA_WEIGHTS_S      9 
#define SQT_DATA_WEIGHTS_D     10 
#define SQT_DATA_INDICES       11
#define SQT_DATA_KNM_CACHE     12
#define SQT_DATA_N_QUAD_POINTS 13

#ifdef SQT_SINGLE_PRECISION

#define SQT_REAL gfloat

#define SQT_FUNCTION_NAME(SQT_func) SQT_func##_f

#define SQRT(SQT_x) sqrtf((SQT_x))
#define CBRT(SQT_x) cbrtf((SQT_x))
#define SIN(SQT_x) sinf((SQT_x))
#define COS(SQT_x) cosf((SQT_x))
#define ACOS(SQT_x) acosf((SQT_x))
#define ATAN(SQT_x) atanf((SQT_x))
#define ATAN2(SQT_y,SQT_x) atan2f((SQT_y),(SQT_x))
#define LOG(SQT_x) logf((SQT_x))
#define EXP(SQT_x) expf(SQT_x)

#else

#define SQT_REAL gdouble

#define SQT_FUNCTION_NAME(SQT_func) SQT_func

#define SQRT(SQT_x) sqrt((SQT_x))
#define CBRT(SQT_x) cbrt((SQT_x))
#define SIN(SQT_x) sin((SQT_x))
#define COS(SQT_x) cos((SQT_x))
#define ACOS(SQT_x) acos((SQT_x))
#define ATAN(SQT_x) atan((SQT_x))
#define ATAN2(SQT_y,SQT_x) atan2((SQT_y),(SQT_x))
#define LOG(SQT_x) log((SQT_x))
#define EXP(SQT_x) exp(SQT_x)

#endif /*SQT_SINGLE_PRECISION*/

#define SIGN(SQT_x) ((SQT_x) < 0 ? -1 : 1)

#define sqt_invert2x2(_Ai, _A)			\
do {						\
  SQT_REAL _det = _A[0]*_A[3] - _A[1]*_A[2] ;	\
						\
  _Ai[0] =  _A[3]/_det ; _Ai[1] = -_A[1]/_det ;	\
  _Ai[2] = -_A[2]/_det ; _Ai[3] =  _A[0]/_det ;	\
 } while (0)
  
#define sqt_vector_scalar(SQT_A,SQT_B)				\
  (((SQT_A)[0])*((SQT_B)[0])+					\
   ((SQT_A)[1])*((SQT_B)[1])+					\
   ((SQT_A)[2])*((SQT_B)[2]))
#define sqt_vector_length(SQT_A)				\
  (SQRT(((SQT_A)[0])*((SQT_A)[0])+				\
	((SQT_A)[1])*((SQT_A)[1]) +				\
	((SQT_A)[2])*((SQT_A)[2])))
#define sqt_vector_length2(SQT_A)				\
  (((SQT_A)[0])*((SQT_A)[0]) +					\
   ((SQT_A)[1])*((SQT_A)[1]) +					\
   ((SQT_A)[2])*((SQT_A)[2]))
#define sqt_vector_cross(SQT_C,SQT_A,SQT_B)				\
  ((SQT_C)[0] = (SQT_A)[1]*(SQT_B)[2] - (SQT_A)[2]*(SQT_B)[1],		\
   (SQT_C)[1] = (SQT_A)[2]*(SQT_B)[0] - (SQT_A)[0]*(SQT_B)[2],		\
   (SQT_C)[2] = (SQT_A)[0]*(SQT_B)[1] - (SQT_A)[1]*(SQT_B)[0])
#define sqt_vector_distance2(SQT_A,SQT_B)		\
  ( ((SQT_A)[0]-(SQT_B)[0])*((SQT_A)[0]-(SQT_B)[0]) +	\
    ((SQT_A)[1]-(SQT_B)[1])*((SQT_A)[1]-(SQT_B)[1]) +	\
    ((SQT_A)[2]-(SQT_B)[2])*((SQT_A)[2]-(SQT_B)[2]) )

#define sqt_vector_diff(SQT_A,SQT_B,SQT_C)	\
  do {						\
  (SQT_A)[0] = (SQT_B)[0] - (SQT_C)[0] ;	\
  (SQT_A)[1] = (SQT_B)[1] - (SQT_C)[1] ;	\
  (SQT_A)[2] = (SQT_B)[2] - (SQT_C)[2] ;	\
  } while (0)

#define sqt_vector_distance(SQT_A,SQT_B)	\
  (SQRT((sqt_vector_distance2(SQT_A,SQT_B))))
#define sqt_vector_diff_scalar(SQT_A,SQT_B,SQT_C)			\
  (((SQT_A)[0]-(SQT_B)[0])*((SQT_C)[0]) +				\
   ((SQT_A)[1]-(SQT_B)[1])*((SQT_C)[1]) +				\
   ((SQT_A)[2]-(SQT_B)[2])*((SQT_C)[2]))

#define sqt_point_interp_jac3(_xe,_xstr,_L0,_L1,_L2,			\
			      _Ls0,_Ls1,_Ls2,_Lt0,_Lt1,_Lt2,_y,_n,_J)	\
  {									\
  (_y)[0] = (_L0)*(_xe)[(_xstr)*0+0] + (_L1)*(_xe)[(_xstr)*1+0] +	\
    (_L2)*(_xe)[(_xstr)*2+0] ;						\
  (_y)[1] = (_L0)*(_xe)[(_xstr)*0+1] + (_L1)*(_xe)[(_xstr)*1+1] +	\
    (_L2)*(_xe)[(_xstr)*2+1] ;						\
  (_y)[2] = (_L0)*(_xe)[(_xstr)*0+2] + (_L1)*(_xe)[(_xstr)*1+2] +	\
    (_L2)*(_xe)[(_xstr)*2+2] ;						\
  (_n)[0] =								\
    ((_Ls0)*(_xe)[(_xstr)*0+1] + (_Ls1)*(_xe)[(_xstr)*1+1] +		\
     (_Ls2)*(_xe)[(_xstr)*2+1])*					\
    ((_Lt0)*(_xe)[(_xstr)*0+2] + (_Lt1)*(_xe)[(_xstr)*1+2] +		\
     (_Lt2)*(_xe)[(_xstr)*2+2]) -					\
    ((_Lt0)*(_xe)[(_xstr)*0+1] + (_Lt1)*(_xe)[(_xstr)*1+1] +		\
     (_Lt2)*(_xe)[(_xstr)*2+1])*					\
    ((_Ls0)*(_xe)[(_xstr)*0+2] + (_Ls1)*(_xe)[(_xstr)*1+2] +		\
     (_Ls2)*(_xe)[(_xstr)*2+2]) ;					\
  (_n)[1] =								\
    ((_Ls0)*(_xe)[(_xstr)*0+2] + (_Ls1)*(_xe)[(_xstr)*1+2] +		\
     (_Ls2)*(_xe)[(_xstr)*2+2])*					\
    ((_Lt0)*(_xe)[(_xstr)*0+0] + (_Lt1)*(_xe)[(_xstr)*1+0] +		\
     (_Lt2)*(_xe)[(_xstr)*2+0]) -					\
    ((_Lt0)*(_xe)[(_xstr)*0+2] + (_Lt1)*(_xe)[(_xstr)*1+2] +		\
     (_Lt2)*(_xe)[(_xstr)*2+2])*					\
    ((_Ls0)*(_xe)[(_xstr)*0+0] + (_Ls1)*(_xe)[(_xstr)*1+0] +		\
     (_Ls2)*(_xe)[(_xstr)*2+0]) ;					\
  (_n)[2] =								\
    ((_Ls0)*(_xe)[(_xstr)*0+0] + (_Ls1)*(_xe)[(_xstr)*1+0] +		\
     (_Ls2)*(_xe)[(_xstr)*2+0])*					\
    ((_Lt0)*(_xe)[(_xstr)*0+1] + (_Lt1)*(_xe)[(_xstr)*1+1] +		\
     (_Lt2)*(_xe)[(_xstr)*2+1]) -					\
    ((_Lt0)*(_xe)[(_xstr)*0+0] + (_Lt1)*(_xe)[(_xstr)*1+0] +		\
     (_Lt2)*(_xe)[(_xstr)*2+0])*					\
    ((_Ls0)*(_xe)[(_xstr)*0+1] + (_Ls1)*(_xe)[(_xstr)*1+1] +		\
     (_Ls2)*(_xe)[(_xstr)*2+1]) ;					\
  (_J) = sqt_vector_length((_n)) ;					\
  (_n)[0] /= (_J) ; (_n)[1] /= (_J) ; (_n)[2] /= (_J) ;			\
  } while (0)

#define sqt_shape_derivatives3(_s,_t,_L,_Ls,_Lt)			\
  {(_L)[0] = 1.0 - (_s) - (_t) ; (_L)[1] = (_s) ; (_L)[2] = (_t) ;	\
  (_Ls)[0] = -1.0 ; (_Ls)[1] =  1.0 ; (_Ls)[2] =  0.0 ;			\
  (_Lt)[0] = -1.0 ; (_Lt)[1] =  0.0 ; (_Lt)[2] =  1.0 ;			\
} while (0)

#define IS_EVEN(SQT_i) (((SQT_i)%2==0)?1:0)

#define minus_one_pow(SQT_n) ((2*((SQT_n)/2) == (SQT_n) ? 1 : -1))

#define yes_if_true(SQT_t)  ((SQT_t) == TRUE ? "yes" : "no") 

#define sqt_cache_position(_xc,_str,_i) (&((_xc)[(_i)*(_str)+0]))
#define sqt_cache_normal(_xc,_str,_i) (&((_xc)[(_i)*(_str)+3]))
#define sqt_cache_weight(_xc,_str,_i) ((_xc)[(_i)*(_str)+6])
#define sqt_cache_coordinate_s(_xc,_str,_i) ((_xc)[(_i)*(_str)+7])
#define sqt_cache_coordinate_t(_xc,_str,_i) ((_xc)[(_i)*(_str)+8])

extern gdouble WANDZURA_7[], WANDZURA_25[], WANDZURA_54[], WANDZURA_85[],
  WANDZURA_126[], WANDZURA_175[], XIAO_GIMBUTAS_453[] ;
extern gfloat WANDZURA_7_F[], WANDZURA_25_F[], WANDZURA_54_F[],
  WANDZURA_85_F[], WANDZURA_126_F[], WANDZURA_175_F[],
  XIAO_GIMBUTAS_453_F[] ;

gint newman_tri(gdouble p[], gdouble x1[], gdouble x2[], gdouble x3[],
		gdouble Iq[], gdouble J[]) ;
gint newman_tri_shape(gdouble p[], gdouble x1[], gdouble x2[], gdouble x3[],
		      gdouble *Imn, gint hmax,
		      gdouble Iq[], gdouble J[]) ;

gint angular_quadrature_select(gint N, gdouble r0, gdouble th0,
			       gdouble **q, gint *nq) ;
gint radial_quadrature_select(gint N, gdouble d, gdouble **q, gint *nq) ;
gint legendre_quadrature_select(gint N, gdouble **q, gint *nq) ;

gint angular_quadrature_select_f(gint N, gfloat r0, gfloat th0,
				 gfloat **q, gint *nq) ;
gint radial_quadrature_select_f(gint N, gfloat d, gfloat **q, gint *nq) ;
gint legendre_quadrature_select_f(gint N, gfloat **q, gint *nq) ;

#endif /*SQT_PRIVATE_H_INCLUDED*/
