function [xstar0,xstar1]=findpoint(a,a1,b,c,d);
%FINDPOINT  
%   [xstar0,xstar1]=findpoint(a,a1,b,c,d) finds points in cx=d that
%   are tangent to x'*Pmid*x=r^2, for some r>0, and to x'*Pup*x=r^2
%
%   Version: 1.0 		    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%


n=size(a,1);

%% Lyapunov function for system A1 (x1star) (middle system)
Pmid=care(a1,zeros(n),eye(n),eye(n),zeros(n),eye(n));
xstar1=d*inv(Pmid)*c'/(c*inv(Pmid)*c');

%% Lyapunov function for system A+Bd (up)
Pup=care(a,zeros(n),eye(n),eye(n),zeros(n),eye(n));
d_in_z=d+c*inv(a)*b*d;
zstar0=d_in_z*inv(Pup)*c'/(c*inv(Pup)*c');
xstar0=zstar0-inv(a)*b*d;


