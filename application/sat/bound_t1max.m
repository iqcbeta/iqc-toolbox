function [tmax]=bound_t1max(a,b,c,tmin,tref);
%BOUND_T1MAX Upper bound the maximum switching times of imapact map 1
%            after a long time. 
%            For Saturation Systems.
%
%   TMAX = BOUND_T1MAX(A,B,C,TMIN,TREF)
%
%   The input TREF is just a reference of time. By default is 20.
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


n=size(a,1);
ainv=inv(a);
F=c;

%%% String variables %%%
A=mat2str(a);
B=mat2str(b);
C=mat2str(c);


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

%time_zero_aprox

%
% Finding exact zeros of Fe^(at)b
%

eq=inline([C '*expm(' A '*t)*' B]);

%fplot(eq,[14 20])
%break

for count=1:num_zeros
	time_zero(count)=fzero(eq,time_zero_aprox(count),optimset('Display','off'));
	count=count+1;
end

%time_zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%	Find t1max  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


findT=0;
T=tmin;

while ~findT
	T=T+dt;

	%%%  calculate  intT

	if num_zeros==0
		sum=abs(F*expm(a*T)*ainv*b);
	else
		count_first=1;
		not_done=1;
		while (count_first<=num_zeros)&(not_done)
			if T<time_zero(count_first)
				not_done=0;
			else
				count_first=count_first+1;
			end
		end
		if count_first==num_zeros+1
			sum=abs(F*expm(a*T)*ainv*b);
		else
			sum=abs(F*(expm(a*time_zero(count_first))-expm(a*T))*ainv*b);
			for count=(count_first+1):num_zeros
				sum=sum+abs(F*(expm(a*time_zero(count))-expm(a*time_zero(count-1)))*ainv*b);
			end
			sum=sum+abs(F*(-expm(a*time_zero(num_zeros)))*ainv*b);
		end
	end
	intT=sum;

	%%%  end calculate  intT

	if (intT+abs(c*expm(a*T)*ainv*b)<=(c*ainv*b+1))
		findT=1;
	end
%T
%intT+abs(c*expm(a*T)*ainv*b)
%pause
end



tmax=T;
