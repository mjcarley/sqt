function J=pjacobi(N,a,b,x)

  J = zeros(N+1, 1) ;

  J(1) = 1.0 ;
  n = 0 ;
  J(n+2) = (2*n+a+b+1)*(2*n+a+b+2)/2/(n+1)/(n+a+b+1)*x + ...
	   (a^2-b^2)*(2*n+a+b+1)/2/(n+1)/(n+a+b+1)/(2*n+a+b) ;

  for n=1:N-1
    J(n+2) = ...
    ((2*n+a+b+1)*(2*n+a+b+2)/2/(n+1)/(n+a+b+1)*x + ...
     (a^2-b^2)*(2*n+a+b+1)/2/(n+1)/(n+a+b+1)/(2*n+a+b))*J(n+1) - ...
    (n+a)*(n+b)*(2*n+a+b+2)/(n+1)/(n+a+b+1)/(2*n+a+b)*J(n) ;
  endfor
  
