function [l1]=l1n(a,b,F,tref);
%L1N  Calculates ||Fe^{At}B||_1
%
%   E = L1N(A,B,F,TREF) calculates the L1 norm of Fe^{At}B.  The input
%   TREF is just a reference of time. By default is 20.
%
%   Version: 1.0		    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%

%%%%%%%%%%  Default variable tref  %%%%%%%%%%%%%%%%
ttt=which('tref');
if isempty(ttt)
        tref=20;
end
clear ttt

%%% String variables %%%
A=mat2str(a);
B=mat2str(b);
Fs=mat2str(F);

n=size(a,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	 l1=||F*expm(a*t)*b||_1         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Finding aprox zeros of Fe^(at)b
%

num_zeros=0;
dt=0.1;
t=dt;
e=1;
antes=F*expm(a*0.01)*b;
while ((e>10^(-7))|(t<tref))
	l1norm=F*expm(a*t)*b;
	if l1norm*antes<0
		num_zeros=num_zeros+1;
		time_zero_aprox(num_zeros)=t;
	end
	e=abs(l1norm-antes);
	antes=l1norm;
	t=t+dt;
end


%
% Finding exact zeros of Fe^(at)b
%

eq=inline([Fs '*expm(' A '*t)*' B]);

for count=1:num_zeros
	time_zero(count)=fzero(eq,time_zero_aprox(count),optimset('Display','off'));
	count=count+1;
end

%
% Finding l1 norm
%
if num_zeros==0
	sum=F*inv(a)*b;
else
	sum=(-1)*F*(expm(a*time_zero(1))-eye(n))*inv(a)*b;
	for count=2:num_zeros
		sum=sum+(-1)^count*F*(expm(a*time_zero(count))-expm(a*time_zero(count-1)))*inv(a)*b;
	end
	sum=sum+(-1)^(num_zeros+1)*F*(-expm(a*time_zero(num_zeros)))*inv(a)*b;
end


l1=abs(sum);


