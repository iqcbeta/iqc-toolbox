function h=factorize(g)
% function h=factorize(g)
% this is a simple spectral factorization routine for
% scalar  polynomials without zeros on the
% imaginary axis: g(w^2)=h(jw)*h(-jw), where
% h is a Hurwitz polynomial

n=length(g)-1;
nn=n:-1:0;
g1=g.*(-1).^nn;
p=zeros(1,2*n+1);
p(1:2:2*n+1)=g1;
z=roots(p);
zp=z(find(real(z)<0));
h=poly(zp);
h=h*sqrt(g(n+1))/h(n+1);

