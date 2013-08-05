function [vperp]=perp(v);
%
%	Given m orthogonal vectors v1...vm, finds n-m unitary vectors
% orthogonal to all vi, that is, vperp'*v=0 (all column vectors of vperp
% are perpendicular to all column vetors of v)
%

[n m]=size(v);

if n>1

% we are given m orthogonal vectors

mmax=m;
for count=1:mmax
	if v(:,count)==0
		m=m-1;
	else
		z(:,count)=v(:,count)/norm(v(:,count));
	end
end

clear count
vperp=[];

count1=m+1;
while count1<=n
	count2=1;
	x=[];
	while count2<=n
		xirand=sign(rand-0.5)*rand;
		x=[x;xirand];
		count2=count2+1;
	end
	y=x;
	for count=1:m
		y=y-(x'*z(:,count))*z(:,count);
	end
	for count3=1:(count1-(m+1))
		y=y-(x'*vperp(:,count3))*vperp(:,count3);
	end
	zvar=y/norm(y);
		%%  to make sure we get always the same vperp for a given v
		%count4=1;
		%while zvar(count4)==0
		%	count4=count4+1;
		%end
		%if sign(zvar(count4))<0
		%	zvar=-zvar;
		%end
	vperp=[vperp zvar];
	count1=count1+1;
end

else
	vperp=v;
end
