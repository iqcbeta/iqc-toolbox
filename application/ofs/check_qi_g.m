function [pas1,pas2]=check_qi_g(a1,iab1,tmax1,xp1,Pi,c,P1,P2,g1,g2,alpha,a,xp0,tmax2);
%%% check if quadratic stab conditions hold for t>tmaxi.  This is a "brute"
%   force way of doing this.  When I have time I will improve this to a more
%   rigorous condition similar to the one when d=0.
%
%   Version: 0.1    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%

pas1=0;
pas2=0;
dterr=0.0000001;

if max(real(eig(a1)))<0
 t1=tmax1;
 vt1=expm(a1*t1)*(xp0+iab1) - iab1 - xp1;
 qt1=c*expm(a1*t1)*Pi;
 pt1=-c*vt1;
 wt1=qt1/pt1;
 Ft1=Pi'*(expm(a1*t1)*Pi+vt1*wt1);
 sigg1=P1-Ft1'*P2*Ft1 - 2*(g1-Ft1'*g2)*wt1 + alpha*wt1'*wt1;
 ea1=min(eig(sigg1+sigg1'));
 err=1;
 t1=tmax1+0.1;
 antes=ea1;
 while err>dterr
	vt1=expm(a1*t1)*(xp0+iab1) - iab1 - xp1;
	qt1=c*expm(a1*t1)*Pi;
	pt1=-c*vt1;
	wt1=qt1/pt1;
	Ft1=Pi'*(expm(a1*t1)*Pi+vt1*wt1);
	sigg1=P1-Ft1'*P2*Ft1 - 2*(g1-Ft1'*g2)*wt1 + alpha*wt1'*wt1;
	ean=min(eig(sigg1+sigg1'));
	if ean<ea1
%	   t1
	   ea1=ean;
%	   pause
	end
	t1=t1+0.1;
	err=abs(antes-ean);
	antes=ean;
 end
% t1
% ea1
 if ea1>0
	pas1=1;
 end
end

if max(real(eig(a)))<0
 t2=tmax2;
 vt2=expm(a*t2)*xp1 - xp0;
 qt2=c*expm(a*t2)*Pi;
 pt2=-c*vt2;
 wt2=qt2/pt2;
 Ft2=Pi'*(expm(a*t2)*Pi+vt2*wt2);
 sigg2=P2-Ft2'*P1*Ft2 - 2*(g2-Ft2'*g1)*wt2 - alpha*wt2'*wt2;
 ea2=min(eig(sigg2+sigg2'));
 err=1;
 t2=tmax2+0.1;
 antes=ea2;
 while err>dterr
	vt2=expm(a*t2)*xp1 - xp0;
	qt2=c*expm(a*t2)*Pi;
	pt2=-c*vt2;
	wt2=qt2/pt2;
	Ft2=Pi'*(expm(a*t2)*Pi+vt2*wt2);
	sigg2=P2-Ft2'*P1*Ft2 - 2*(g2-Ft2'*g1)*wt2 - alpha*wt2'*wt2;
	ean=min(eig(sigg2+sigg2'));
	if ean<ea2
%	   t2
	   ea2=ean;
%	   pause
	end
	t2=t2+0.1;
	err=abs(antes-ean);
	antes=ean;
 end
% t2
% ea2
 if ea2>0
	pas2=1;
 end
end

