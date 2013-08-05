function [tmin]=bound_tmin(a,b,c,d,tref);
%BOUND_TMIN Lower bound on the minimum switching time after a long time
%           for hysterisis relay.
%
%   TMIN = BOUND_TMIN(A,B,C,D,TREF)
%
%   The input TREF is just a reference of time. By default is 20.
%
%   Version: 1.1    Last modified: Aug 9, 1999
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%

%%%%%%%%%%  Default variable d  %%%%%%%%%%%%%%%%
ttt=which('tref');
if isempty(ttt)
        tref=20;
end
clear ttt

n=size(a,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	Find tmin1     %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Finding kdd = kdd1 + kdd2
%

kdd1=l1n(a,b,c*a^2,tref);

%       kdd2=max_{t} abs(c*expm(a*t)*a*b)

kdd2=abs(c*a*b);
antes=kdd2;
dt=0.1;
t=dt;
e=1;
while ((e>10^(-4))|(t<tref))
	kdd2try=abs(c*expm(a*t)*a*b);
	if kdd2try>kdd2
		kdd2=kdd2try;
	end
	e=abs(kdd2try-antes);
	antes=kdd2try;
	t=t+dt;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kdd = kdd1 + kdd2;
kd=-2*c*b;

tmin1=(kd+sqrt(kd^2+4*kdd*d))/kdd;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	Find tmin2     %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Finding kdl = kdl1 + kdl2
%

kdl1=l1n(a,b,c*a,tref);

%       kdl2=max_{t} abs(c*expm(a*t)*b)

kdl2=abs(c*b);
antes=kdl2;
dt=0.1;
t=dt;
e=1;
while ((e>10^(-4))|(t<tref))
	kdl2try=abs(c*expm(a*t)*b);
	if kdl2try>kdl2
		kdl2=kdl2try;
	end
	e=abs(kdl2try-antes);
	antes=kdl2try;
	t=t+dt;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kdl = kdl1 + kdl2;

tmin2=2*d/kdl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmin = max([tmin1 tmin2]);

