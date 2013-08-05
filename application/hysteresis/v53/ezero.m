function [all_zeros]=ezero(g,tmin,tmax,dt);
%EZERO   Find all zeros of g between tmin and tmax expect at tmin.



%%%%%%%%%%  Default variable d  %%%%%%%%%%%%%%%%
ttt=which('dt');
if isempty(dt)
        dt=0.01;
end
clear ttt


%
% Finding aprox zeros of g
%

num_zeros=0;
t=tmin+dt;
antes=g(tmin+0.000001);
while t<tmax
	gfun=g(t);
	if gfun*antes<0
		num_zeros=num_zeros+1;
		time_zero_aprox(num_zeros)=t;
	end
	antes=gfun;
	t=t+dt;
end

%time_zero_aprox
%num_zeros

%
% Finding exact zeros of g
%

all_zeros=[];
for count=1:num_zeros
	all_zeros(count)=fzero(g,time_zero_aprox(count),optimset('Display','off'));
end

