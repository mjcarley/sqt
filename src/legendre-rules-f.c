/* automatically-generated file, do not edit */
/* Thu Dec  3 12:59:15 GMT 2020 */

#include <glib.h>

gfloat legendre_rule4_f[] = {
-7.7459666924148340e-01, 5.5555555555555536e-01,
0.0000000000000000e+00, 8.8888888888888884e-01,
7.7459666924148340e-01, 5.5555555555555536e-01,
0.0} ;
gfloat legendre_rule8_f[] = {
-9.0617984593866396e-01, 2.3692688505618917e-01,
-5.3846931010568311e-01, 4.7862867049936664e-01,
0.0000000000000000e+00, 5.6888888888888889e-01,
5.3846931010568311e-01, 4.7862867049936664e-01,
9.0617984593866396e-01, 2.3692688505618917e-01,
0.0} ;
gfloat legendre_rule12_f[] = {
-9.4910791234275860e-01, 1.2948496616886979e-01,
-7.4153118559939446e-01, 2.7970539148927676e-01,
-4.0584515137739718e-01, 3.8183005050511903e-01,
0.0000000000000000e+00, 4.1795918367346940e-01,
4.0584515137739718e-01, 3.8183005050511903e-01,
7.4153118559939446e-01, 2.7970539148927676e-01,
9.4910791234275860e-01, 1.2948496616886979e-01,
0.0} ;
gfloat legendre_rule16_f[] = {
-9.6816023950762609e-01, 8.1274388361574690e-02,
-8.3603110732663577e-01, 1.8064816069485815e-01,
-6.1337143270059036e-01, 2.6061069640293605e-01,
-3.2425342340380892e-01, 3.1234707704000286e-01,
0.0000000000000000e+00, 3.3023935500125978e-01,
3.2425342340380892e-01, 3.1234707704000286e-01,
6.1337143270059036e-01, 2.6061069640293605e-01,
8.3603110732663577e-01, 1.8064816069485815e-01,
9.6816023950762609e-01, 8.1274388361574690e-02,
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
