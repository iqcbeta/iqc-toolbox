function not_origin_c(G,d,nu1,nu2,tmax1,tmax2,npoi1,npoi2);
%
% Stability of on_off systems when d ~= 0 (equilibrium point does not
%  belong to the switching surface).  Less conservative method.  Here, 
%  we use the S-Procedure.
%
%   Version: 0.1    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%


%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu



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
else
    disp(' ')
    disp(' d=0.  Use stability condition for d=0 in directory "with_origin"')
    disp(' ')
    break
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

[xu,xd]=findpoint(a,a1,b1,c,d,pa1);
xp0=xu;
xp1=xd;


%%%%%%%%%   Boundary of S_ai   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (same for S_a1 and S_a2, but different when seen from xp0 or xp1)

q0=c*a*Pi;
p01=-c*a*xp0;
l01=q0'*p01/(q0*q0');

p02=-c*a*xp1;
l02=q0'*p02/(q0*q0');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  LMI STABILITY CONDITIONS  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmin1=0.001;
tmin2=0.001;

abst_init_lmi;

lmitbx_options([0 0 0 0 1]);

P1=symmetric(n-1);
P2=symmetric(n-1);
G1=diagonal(n-1);
G2=diagonal(n-1);
alpha=symmetric;
tau1=diagonal(nu1);
tau2=diagonal(nu2);

P1>0;
P2>0;
g1=G1*ones(n-1,1);
g2=G2*ones(n-1,1);
tau1>0;
tau2>0;


dt1=(tmax1-tmin1)/(nu1-1);
timee1=tmin1:dt1:tmax1;
count=1;
for t1=timee1
	vt1=expm(a1*t1)*(xp0+iab1) - iab1 - xp1;
	qt1=c*expm(a1*t1)*Pi;
	pt1=-c*vt1;
	wt1=qt1/pt1;
	lttemp=qt1'*pt1/(qt1*qt1');
	sttemp=p01*(l01/norm(l01)*norm(lttemp)-lttemp/norm(lttemp)*norm(l01));
	sttemp=sttemp/norm(sttemp);
	betat=lttemp*sttemp'+sttemp*lttemp';
	lt1(:,count)=lttemp;
	st1(:,count)=sttemp;
	Ft1=Pi'*(expm(a1*t1)*Pi+vt1*wt1);
	P1-Ft1'*P2*Ft1 - 2*(g1-Ft1'*g2)*wt1 + alpha*wt1'*wt1 - tau1(count,count)*betat >0;
	count=count+1;
end

clear t1

dt2=(tmax2-tmin2)/(nu2-1);
timee2=tmin2:dt2:tmax2;
count=1;
for t2=timee2
	vt2=expm(a*t2)*xp1 - xp0;
	qt2=c*expm(a*t2)*Pi;
	pt2=-c*vt2;
	wt2=qt2/pt2;
	lttemp=qt2'*pt2/(qt2*qt2');
	sttemp=-p02*(l02/norm(l02)*norm(lttemp)-lttemp/norm(lttemp)*norm(l02));
	sttemp=sttemp/norm(sttemp);
	betat=lttemp*sttemp'+sttemp*lttemp';
	lt2(:,count)=lttemp;
	st2(:,count)=sttemp;
	Ft2=Pi'*(expm(a*t2)*Pi+vt2*wt2);
	P2-Ft2'*P1*Ft2 - 2*(g2-Ft2'*g1)*wt2 - alpha*wt2'*wt2 - tau2(count,count)*betat >0;
	count=count+1;
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
tau1=value_iqc(tau1);
tau2=value_iqc(tau2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plotting min(eig(Qi-Ft'*Q*Ft)) for tmin<t<tmax
%

count=1;
for t1=timee1	
	vt1=expm(a1*t1)*(xp0+iab1) - iab1 - xp1;
	qt1=c*expm(a1*t1)*Pi;
	pt1=-c*vt1;
	wt1=qt1/pt1;
	betat=lt1(:,count)*st1(:,count)'+st1(:,count)*lt1(:,count)';
	Ft1=Pi'*(expm(a1*t1)*Pi+vt1*wt1);
	sigg1=P1-Ft1'*P2*Ft1 - 2*(g1-Ft1'*g2)*wt1 + alpha*wt1'*wt1 -tau1(count,count)*betat;
	Rstab1(count)=min(eig(sigg1+sigg1'))/2;
	count=count+1;
end

count=1;
for t2=timee2
	vt2=expm(a*t2)*xp1 - xp0;
	qt2=c*expm(a*t2)*Pi;
	pt2=-c*vt2;
	wt2=qt2/pt2;
	betat=lt2(:,count)*st2(:,count)'+st2(:,count)*lt2(:,count)';
	Ft2=Pi'*(expm(a*t2)*Pi+vt2*wt2);
	sigg2=P2-Ft2'*P1*Ft2 - 2*(g2-Ft2'*g1)*wt2 - alpha*wt2'*wt2 -tau2(count,count)*betat;
	Rstab2(count)=min(eig(sigg2+sigg2'))/2;
	count=count+1;
end

figure(2)
clf
plot(timee1,Rstab1,'b')
hold on
plot(timee2,Rstab2,'r')
grid on
title('Min(eig(R_i(t))); Want them always positive')

%min(Rstab1);
%min(Rstab2);

if min(Rstab1)>0 & min(Rstab2)>0 & n>3
    check_qi_g
end   

pas1=0;
pas2=0;
if min(Rstab1)>-1e-12 & min(Rstab2)>-1e-12
    [pas1,pas2]=check_qi_g(a1,iab1,tmax1,xp1,Pi,c,P1,P2,g1,g2,alpha,a,xp0,tmax2);
end   



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
        disp(' Increase tmax1 and tmax2 and run "ofs" again')
        disp(' ')
        disp(' ')
end

