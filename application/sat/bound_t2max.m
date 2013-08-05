function [tmax]=bound_t2max(a,b,c,tmin,dt,tref);
%BOUND_T2MAX Upper bound the maximum switching times of imapact maps 2a
%	     and 2b after a long time. 
%            For Saturation Systems.
%
%   TMAX = BOUND_T2MAX(A,B,C,TMIN,dt,TREF)
%
%   dt is the precision
%   The input TREF is just a reference of time. By default is 20.
%
%   Version: 0.1		    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%
%   Note: need to improve efficiency of this function!!
%

%%%%%%%%%%  Default variable tref  %%%%%%%%%%%%%%%%
ttt=which('tref');
if isempty(ttt)
        tref=20;
end
clear ttt

%%%%%%%%%%  Default variable dt  %%%%%%%%%%%%%%%%
ttt=which('dt');
if isempty(ttt)
        dt=0.1;
end
clear ttt


n=size(a,1);
a1=a+b*c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%	Find t1max  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


findT=0;
T=tmin;

while ~findT
	if T<1
	  T=T+dt;
	else
	  T=T*(1+dt);
	end
	F=c*expm(a1*T);
	T1n=l1n(a,b,F,tref);
	if T1n<=1
		findT=1;
	end
%T
%T1n
%pause
end



tmax=T;
