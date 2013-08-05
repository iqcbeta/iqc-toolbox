function [xstar0,xstar1]=findpoint(a,a1,b1,c,d,pa1);
%FINDPOINT  
%   [xstar0,xstar1]=findpoint(a,a1,b1,c,d) finds points in cx=d that
%   are tangent to x'*Pdown*x=r^2, for some r>0, and to x'*Pup*x=r^2
%
%   Version: 0.1    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%


n=size(a,1);
iab1=inv(a1)*b1;

%% Lyapunov function of each system A (down)
Pdown=care(a,zeros(n),eye(n),eye(n),zeros(n),eye(n));
xstar1=d*inv(Pdown)*c'/(c*inv(Pdown)*c');

%% Need to make sure A1 is stable before finding a Lyapunov function
if max(real(pa1))<0
  Pup=care(a1,zeros(n),eye(n),eye(n),zeros(n),eye(n));
  d_in_z=d+c*inv(a1)*b1;
  zstar0=d_in_z*inv(Pup)*c'/(c*inv(Pup)*c');
  xstar0=zstar0-inv(a1)*b1;
else
  disp(' A1 has complex unstable poles.  Continuing... ')
  disp(' ')
     %%%% going to find xstar0 along a stable mode
  [v,pa1]=eig(a1);
  notfound=1;
  m=1;
  while notfound
     if isreal(pa1(m,m))
	% real stable eigenvalue
	xstar0=-iab1+(d+c*iab1)/(c*v(:,m))*v(:,m);
	notfound=0;
     elseif real(pa1(m,m))<0
	% complex stable eigenvalue       % see Feb 25, 2k
	va=v(:,m)+conj(v(:,m));
	va=va/norm(va);
	vb=i*(v(:,m)-conj(v(:,m)));
	vb=vb/norm(vb);
	vc=-(va'*vb)*va+vb;
	vc=vc/norm(vc);
	V=[va vc];  %% orthogonal stable basis
	av=V'*a1*V;
	Pv=care(av,zeros(2),eye(2),eye(2),zeros(2),eye(2));
	alpha=(d+c*iab1)/(c*V*inv(Pv)*V'*c')*inv(Pv)*V'*c';
	xstar0=-iab1+V*alpha;
	notfound=0;
     else
       m=m+2;
     end
  end
end


