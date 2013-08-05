function [tmin]=bound_t2bmin(a,b,c,tref);
%BOUND_TMIN Lower bound on the minimum switching time of impact map 2b
%	     after a long time
%           For SAT.
%
%   TMIN = BOUND_2BTMIN(A,B,C,D,TREF)
%
%   The input TREF is just a reference of time. By default is 20.
%
%   Version: 1.0		    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%

%%%%%%%%%%  Default variable d  %%%%%%%%%%%%%%%%
ttt=which('tref');
if isempty(ttt)
        tref=20;
end
clear ttt

n=size(a,1);
a1=a+b*c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	Find tmin1     %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kdd=l1n(a,b,c*a1^2,tref);
tmin1=2/sqrt(kdd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%	Find tmin2     %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kdl=l1n(a,b,c*a1,tref);
tmin2=2/kdl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmin = max([tmin1 tmin2]);

