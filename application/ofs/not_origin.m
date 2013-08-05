function not_origin(G,d,nu1,nu2,tmax1,tmax2,npoi1,npoi2);
%
% Stability of on_off systems when d ~= 0 (equilibrium point does not
%  belong to the switching surface)
%
%   Version: 0.1    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%
%%% This version does not include the S-procedure. 
%%% Use not_origin_c to include conic a relation.  This version includes
%%% conditions at 0 (special structure for P2 and g2). (This conditions
%%% restricts this to work with n>2 only!)


%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu



%clear
%s=tf([1 0],1);
%G=-(2*s^2+2*s+12)/(s^3+2*s^2+2*s+3);
%d=1;
%nu=40;

[at,bt,c]=ssdata(G);
dt=d;
a1t=at+bt*c;
b1t=-bt*dt;

%make sure system has one eq point
if d>0
  if -c*inv(a1t)*b1t>=dt
    disp(' ')
    disp(' System has two equilibrium points (d>0)')
    disp(' ')
    break
  else
    a=at; d=dt; a1=a1t; b1=b1t;
  end 
elseif d<0
  if -c*inv(a1t)*b1t<=dt
    disp(' ')
    disp(' System has no equilibrium points (d<0)')
    disp(' ')
    break
  else
    a=a1t; d=-(dt+c*inv(a1t)*b1t); a1=at; b1=at*inv(a1t)*b1t;
  end
end
  

%% check if origin is locally stable and if there exists a 
%% real unstable pole of a1
stop=0;
pa=eig(a);
for cc=1:1:size(pa,1)
   if real(pa(cc))>=0
      disp(' ')
      disp(' The origin is locally unstable (A has unstable poles)')
      disp(' ')
      stop=1;
      break
   end
end

pa1=eig(a1);
for cc=1:1:size(pa1,1)
   if imag(pa1(cc))==0 & pa1(cc)>0
      disp(' ')
      disp(' The system is unstable (A+BC has unstable real poles)')
      disp(' ')
      stop=1;
      break
   end
end
if stop
   break
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% End of checks.  Begining of the program %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=size(a,1);
Pi=perp(c');
iab1=inv(a1)*b1;


%%%%%%%%%%%%%%%%  Get xp0 and xp1  %%%%%%%%%%%%%

[xstar0,xstar1]=findpoint(a,a1,b1,c,d,pa1);
xp0=xstar0;
xp1=xstar1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  LMI STABILITY CONDITIONS  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmin1=0.01;
tmin2=tmin1;

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
w0=-(c*a*Pi)/(c*a*xp0);
l=perp(w0');
z=w0'/(w0*w0');
F=[l z/norm(z)];
P1=F*[zeros(n-2) lam1; lam1' lam2]*F' + P2;

%%   define g1
v0=Pi'*(xp0-xp1);
g1 = F*[lam1;lam2]*norm(z) - P2*v0 + g2 - k*z;

%%   define alpha

alpha= lam2*z'*z + - 2*k*z'*z + v0'*P2*v0 - 2*v0'*g2;

P1>0;
P2>0;

dt1=(tmax1-tmin1)/(nu1-1);
timee1=tmin1:dt1:tmax1;
for t1=timee1
	vt1=expm(a1*t1)*(xp0+iab1) - iab1 - xp1;
	qt1=c*expm(a1*t1)*Pi;
	pt1=-c*vt1;
	wt1=qt1/pt1;
	Ft1=Pi'*(expm(a1*t1)*Pi+vt1*wt1);
	P1-Ft1'*P2*Ft1 - 2*(g1-Ft1'*g2)*wt1 + wt1'*alpha*wt1 >0;
end

clear t1

dt2=(tmax2-tmin2)/(nu2-1);
timee2=tmin2:dt2:tmax2;
for t2=timee2
	vt2=expm(a*t2)*xp1 - xp0;
	qt2=c*expm(a*t2)*Pi;
	pt2=-c*vt2;
	wt2=qt2/pt2;
	Ft2=Pi'*(expm(a*t2)*Pi+vt2*wt2);
	P2-Ft2'*P1*Ft2 - 2*(g2-Ft2'*g1)*wt2 - wt2'*alpha*wt2 >0;
end


lmi_mincx_tbx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  End of Main LMI's  definition  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P1=value_iqc(P1);
P2=value_iqc(P2);
g1=value_iqc(g1);
g2=value_iqc(g2);
alpha=value_iqc(alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plotting min(eig(Qi-Ft'*Q*Ft)) for tmin<t<tmax
%

ex1=5;
ex2=ex1;
clear timee1 dt1;  dt1=(tmax1+ex1)/npoi1; timee1=0:dt1:(tmax1+ex1);
%timee1=[0 timee1];
count=1;
for t1=timee1	
	if t1==0
	  vt1=xp0-xp1;
	  wt1=-(c*a1*Pi)/(c*a1*(xp0+iab1));
	  Ft1=Pi'*(Pi+vt1*wt1);
	else
	  vt1=expm(a1*t1)*(xp0+iab1) - iab1 - xp1;
	  qt1=c*expm(a1*t1)*Pi;
	  pt1=-c*vt1;
	  wt1=qt1/pt1;
	  Ft1=Pi'*(expm(a1*t1)*Pi+vt1*wt1);
	end
	sigg1=P1-Ft1'*P2*Ft1 - 2*(g1-Ft1'*g2)*wt1 + alpha*wt1'*wt1;
	Rstab1(count)=min(eig(sigg1+sigg1'))/2;
	count=count+1;
end

clear timee2 dt2;  dt2=(tmax2+ex2)/npoi2; timee2=0:dt2:(tmax2+ex2);
%timee2=[0 timee2];
count=1;
for t2=timee2
	if t2==0
	  vt2=xp1 - xp0;
	  wt2=-(c*a*Pi)/(c*a*xp1);
	  Ft2=Pi'*(Pi+vt2*wt2);
	else
	  vt2=expm(a*t2)*xp1 - xp0;
	  qt2=c*expm(a*t2)*Pi;
	  pt2=-c*vt2;
	  wt2=qt2/pt2;
	  Ft2=Pi'*(expm(a*t2)*Pi+vt2*wt2);
	end
	sigg2=P2-Ft2'*P1*Ft2 - 2*(g2-Ft2'*g1)*wt2 - alpha*wt2'*wt2;
	Rstab2(count)=1*min(eig(sigg2+sigg2'))/2;
	count=count+1;
end

figure(1)
clf
plot(timee1,Rstab1,'b')
hold on
plot(timee2,Rstab2*1,'r')
grid on
title('Min(eig(R_i(t))); Want them always positive')

%min(Rstab1)
%min(Rstab2)

pas1=0;
pas2=0;
if min(Rstab1)>-1e-12 & min(Rstab2)>-1e-12
    [pas1,pas2]=check_qi_g(a1,iab1,tmax1,xp1,Pi,c,P1,P2,g1,g2,alpha,a,xp0,tmax2);
end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  FINAL CONCLUSIONS  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (pas1 & pas2)
        disp(' ')
        disp(' ')
        disp('   THE EQUILIBRIUN POINT IS GLOBALLY ASYMPTOTICALLY STABLE')
        disp(' ')
        disp(' ')
else
        disp(' ')
        disp(' ')
        disp('   The equilibrium point could not be proven to be globally stable')
        disp(' Increase tmax1 and tmax2 and run "ofs" again, or try the less')
        disp('conservative method using the S-procedure (use 1 in the last entry of ofs)')
        disp(' ')
        disp(' ')
end
