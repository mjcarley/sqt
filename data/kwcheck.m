load ktest.dat

n = ktest(:,1) ;
m = ktest(:,2) ;
u = ktest(:,3) ;
v = ktest(:,4) ;

x = (2*u + v - 1)./(1-v) ;

Kc = ktest(:,5) ;

Kr = 0*Kc ;
for i=1:length(n)
  P = legendre(m(i), x(i)) ;
  P = P(1) ;
  J = pjacobi(n(i)-m(i), 0, 2*m(i)+1, 1-2*v(i)) ;
  J = J(n(i)-m(i)+1) ;
  Kr(i) = (1-v(i))^m(i)*P*J*sqrt(2*(1+2*m(i))*(n(i)+1)) ;
endfor
