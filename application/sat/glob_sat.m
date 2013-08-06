function glob_sat(G,d,nu1,nu2a,nu2b,npoi1,npoi2a,npoi2b);
%
% Stability of saturation systems
%
%   Version: 1.0		    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%
%%% Global stability of a saturations system.   See README
%%% file for more information.  This version does not include the S-procedure. 
%%% A version with glob_sat_c will be available in future versions.
%%% This version includes conditions at 0 (special structure for P2 
%%% and g2). (This conditions restricts this to work with n>2 only!).

%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu



[a,b,c]=ssdata(G);
a1=a+b*c;

if d==0
   disp('  System is linear.  ')
   disp(' ')
%    break
elseif d<0
   d=-d;
end

% A cannot have eigenvalues with positive real part
pa=eig(a);
if max(real(pa))>0
      disp(' ')
      disp(' The system is unstable (A has unstable poles)')
      disp(' ')
%       break
end

%make sure a and a+bc are stale
pa1=eig(a1);
if max(real(pa1))>=0
      disp(' ')
      disp(' The origin is locally unstable (A+BC has unstable poles)')
      disp(' ')
%       break
end


iabd=inv(a)*b*d;

%%%make sure system has only one eq point
if -c*iabd>=d
    disp(' ')
    disp(' System has 3 equilibrium points')
    disp(' ')
%     break
end
  
%break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% End of checks.  Begining of the program %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=size(a,1);
Pi=perp(c');

%%%%%%%%%%%%%%%%  Get xp0 and xp1  %%%%%%%%%%%%%

[xstar0,xstar1]=findpoint(a,a1,b,c,d);
xp0=xstar0;
xp1=xstar1;

%%%%%%%%%%%%%%%%  Bounds on switching times  %%%%%%%%%%%%%
tmin1=0.01; 
tmin2a=tmin1;

disp(' Finding bounds on switching times')
tmin2b=bound_t2bmin(a,b,c);
tmax1=bound_t1max(a,b,c,0);
tmax2a=bound_t2max(a,b,c,tmin2b,0.5);
tmax2b=tmax2a;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  LMI STABILITY CONDITIONS  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' Solving LMIs')

abst_init_lmi;
lmitbx_options([0 0 0 0 1]);

P2=symmetric(n-1);
G2=diagonal(n-1);
g2=G2*ones(n-1,1);

LAM1=diagonal(n-2);
lam1=LAM1*ones(n-2,1);
lam2=symmetric;

k=symmetric;

%%   define P1
w0=-(c*a*Pi)/(c*(a*xp0+b*d));
l=perp(w0');
z=w0'/(w0*w0');
F=[l z/norm(z)];
P1=F*[zeros(n-2) lam1; lam1' lam2]*F' + P2;

%%   define g1
v0=Pi'*(xp0-xp1);
%g1=P1*z-P2*z-P2*v0+g2-k*z;
g1 = F*[lam1;lam2]*norm(z) - P2*v0 + g2 - k*z;

%%   define alpha

alpha= lam2*z'*z + - 2*k*z'*z + v0'*P2*v0 - 2*v0'*g2;

P1>0;
P2>0;

dt1=(tmax1-tmin1)/(nu1-1);
timee1=tmin1:dt1:tmax1;
for t1=timee1
	vt1=expm(a*t1)*(xp0+iabd) - iabd - xp1;
	qt1=c*expm(a*t1)*Pi;
	pt1=-c*vt1;
	wt1=qt1/pt1;
	Ft1=Pi'*(expm(a*t1)*Pi+vt1*wt1);
	P1-Ft1'*P2*Ft1 - 2*(g1-Ft1'*g2)*wt1 + wt1'*alpha*wt1 >0;
end
clear t1


dt2a=(tmax2a-tmin2a)/(nu2a-1);
timee2a=tmin2a:dt2a:tmax2a;
for t2a=timee2a
	vt2a=expm(a1*t2a)*xp1 - xp0;
	qt2a=c*expm(a1*t2a)*Pi;
	pt2a=-c*vt2a;
	wt2a=qt2a/pt2a;
	Ft2a=Pi'*(expm(a1*t2a)*Pi+vt2a*wt2a);
	P2-Ft2a'*P1*Ft2a - 2*(g2-Ft2a'*g1)*wt2a - wt2a'*alpha*wt2a >0;
end
clear t2a vt2a qt2a pt2a wt2a Ft2a


dt2b=(tmax2b-tmin2b)/(nu2b-1);
timee2b=tmin2b:dt2b:tmax2b;
for t2b=timee2b
	vt2b=expm(a1*t2b)*xp1 + xp0;
	qt2b=c*expm(a1*t2b)*Pi;
	pt2b=-c*vt2b;
	wt2b=qt2b/pt2b;
	Ft2b=Pi'*(expm(a1*t2b)*Pi+vt2b*wt2b);
	P2-Ft2b'*P1*Ft2b - 2*(g2+Ft2b'*g1)*wt2b - wt2b'*alpha*wt2b >0;
end
clear t2b vt2b qt2b pt2b wt2b Ft2b


lmi_mincx_tbx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  End of Main LMI's  definition  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' Plotting')

P1=value_iqc(P1);
P2=value_iqc(P2);
g1=value_iqc(g1);
g2=value_iqc(g2);
alpha=value_iqc(alpha);


% Plotting min(eig(Qi-Ft'*Q*Ft)) for tmin<t<tmax  for g_sat
%


clear Rstab1 Rstab2a Rstab2b

ex=0;

clear timee1 dt1;  dt1=(tmax1+ex)/npoi1; timee1=0:dt1:(tmax1+ex);
%timee1=[0 timee1];
count=1;
for t1=timee1	
	if t1==0
	  vt1=xp0-xp1;
	  wt1=-(c*a*Pi)/(c*(a*xp0+b*d));
	  Ft1=Pi'*(Pi+vt1*wt1);
	else
	  vt1=expm(a*t1)*(xp0+iabd) - iabd - xp1;
	  qt1=c*expm(a*t1)*Pi;
	  pt1=-c*vt1;
	  wt1=qt1/pt1;
	  Ft1=Pi'*(expm(a*t1)*Pi+vt1*wt1);
	end
	sigg1=P1-Ft1'*P2*Ft1 - 2*(g1-Ft1'*g2)*wt1 + wt1'*alpha*wt1;
	Rstab1(count)=min(eig(sigg1+sigg1'))/2;
%	Rstabm1(count)=max(eig(sigg1+sigg1'))/2;
	count=count+1;
end

%ex=5;
clear timee2a dt2a;  dt2a=(tmax2a+ex)/npoi2a; timee2a=0:dt2a:(tmax2a+ex);
%timee2a=[0 timee2a];
count=1;
for t2a=timee2a
	if t2a==0
	  vt2a=xp1 - xp0;
	  wt2a=-(c*a1*Pi)/(c*a1*xp1);
	  Ft2a=Pi'*(Pi+vt2a*wt2a);
	else
	  vt2a=expm(a1*t2a)*xp1 - xp0;
	  qt2a=c*expm(a1*t2a)*Pi;
	  pt2a=-c*vt2a;
	  wt2a=qt2a/pt2a;
	  Ft2a=Pi'*(expm(a1*t2a)*Pi+vt2a*wt2a);
	end
	sigg2=P2-Ft2a'*P1*Ft2a - 2*(g2-Ft2a'*g1)*wt2a - wt2a'*alpha*wt2a;
	Rstab2a(count)=1*min(eig(sigg2+sigg2'))/2;
%	Rstabm2a(count)=1*max(eig(sigg2+sigg2'))/2;
	count=count+1;
end
clear t2a vt2a qt2a pt2a wt2a Ft2a


%ex=5;
clear timee2b dt2b;  dt2b=(tmax2b+ex)/npoi2b; timee2b=tmin2b:dt2b:(tmax2b+ex);
%timee2b=[0 timee2b];
count=1;
for t2b=timee2b
%	if t2b==0
%	  vt2b=xp1 + xp0;
%	  wt2b=(c*a1*Pi)/(c*a1*xp1);
%	  Ft2b=Pi'*(-Pi - vt2b*wt2b);
%	else
	  vt2b=expm(a1*t2b)*xp1 + xp0;
	  qt2b=c*expm(a1*t2b)*Pi;
	  pt2b=-c*vt2b;
	  wt2b=qt2b/pt2b;
	  Ft2b=Pi'*(expm(a1*t2b)*Pi+vt2b*wt2b);
%	end
	sigg2=P2-Ft2b'*P1*Ft2b - 2*(g2+Ft2b'*g1)*wt2b - wt2b'*alpha*wt2b;
	Rstab2b(count)=1*min(eig(sigg2+sigg2'))/2;
%	Rstabm2b(count)=1*max(eig(sigg2+sigg2'))/2;
	count=count+1;
end
clear t2b vt2b qt2b pt2b wt2b Ft2b


figure(1)
clf
plot(timee1,Rstab1,'b')
hold on
%plot(timee1,Rstabm1,'b--')
plot(timee2a,Rstab2a,'r')
%plot(timee2a,Rstabm2a,'r--')
plot(timee2b,Rstab2b,'k')
%plot(timee2b,Rstabm2b,'k--')
grid on
title('min(eig(R_i(t))).  Want them to be always positive')

%min(Rstab1)
%min(Rstab2a)
%min(Rstab2b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  FINAL CONCLUSIONS  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if min(Rstab1)>-1e-12 & min(Rstab2a)>-1e-12 & min(Rstab2b)>-1e-12
        disp(' ')
        disp(' ')
        disp('   THE EQUILIBRIUN POINT IS GLOBALLY ASYMPTOTICALLY STABLE')
        disp(' ')
        disp(' ')
else
        disp(' ')
        disp(' ')
        disp('   The equilibrium point could not be proven to be globally stable')
        disp(' Try increasing nu1, nu2a, or nu2b')
        disp(' ')
        disp(' ')
end
