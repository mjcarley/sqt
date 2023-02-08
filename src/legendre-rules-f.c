/* automatically-generated file, do not edit */
/* Fri Apr 16 15:02:53 BST 2021 */

#include <glib.h>

gfloat legendre_rule4_f[] = {
    .11270166537925830000, .03130601816090507231,
    .50000000000000000, .22222222222222221000,
    .88729833462074170000, .24647175961687260768,
0.0} ;
gfloat legendre_rule8_f[] = {
    .04691007703066802000, .00555712921431103087,
    .23076534494715844500, .05522545512469309087,
    .50000000000000000, .14222222222222222250,
    .76923465505284155500, .18408888012499022912,
    .95308992296933198000, .11290631331378355412,
0.0} ;
gfloat legendre_rule12_f[] = {
    .02544604382862070000, .00164744006214026462,
    .12923440720030277000, .01807378022992264681,
    .29707742431130141000, .05671654396435744667,
    .50000000000000000, .10448979591836735000,
    .70292257568869859000, .13419848128820206832,
    .87076559279969723000, .12177891551471573318,
    .97455395617137930000, .06309504302229463037,
0.0} ;
gfloat legendre_rule16_f[] = {
    .01591988024618695500, .00064693926489917993,
    .08198444633668211500, .00740516971815396264,
    .19331428364970482000, .02518988504329214359,
    .33787328829809554000, .05276686700490217221,
    .50000000000000000, .08255983875031494500,
    .66212671170190446000, .10340667151509925778,
    .80668571635029518000, .10511546315817588140,
    .91801555366331788500, .08291891062927511235,
    .98408011975381304500, .03999025491588816506,
0.0} ;


/*basic stub file for Legendre quadrature selection*/

gint legendre_quadrature_select_f(gint N, gfloat **q, gint *nq)

{
  /*snap N to one of the available values*/
  N = 4*(N/4) ;
  if ( N < 4 ) N = 4 ;
  if ( N > 16 ) N = 16 ;

  switch ( N ) {
  case  4:
    *q = legendre_rule4_f ;
    *nq = 3 ;
    break ;
  case  8:
    *q = legendre_rule8_f ;
    *nq = 5 ;
    break ;
  case  12:
    *q = legendre_rule12_f ;
    *nq = 7 ;
    break ;
  case  16:
    *q = legendre_rule16_f ;
    *nq = 9 ;
    break ;
  default: g_assert_not_reached() ; break ; 
  }

  return 0 ;
}
