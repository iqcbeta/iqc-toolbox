function origin(G,nu1,nu2,tmax1,tmax2,npoi1,npoi2);
%
% Stability of on_off systems when d = 0 (equilibrium point belongs
% to the switching surface)
%
%   Version: 0.1    Last modified: March 3, 2001
%   Written by Jorge Goncalves, jmg@mit.edu, All rights reserved
%
%%% This version includes conditions at 0 (special structure for P2 and g2).
%%% (This conditions restricts this to work with n>2 only!)

%
% change call of "value" to "_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu


[a,b,c]=ssdata(G);
n=size(a,1);
Pi=perp(c');
a1=a+b*c;

%% check if there exists a real unstable pole of a and a1
stop=0;
pa=eig(a);
pa1=eig(a1);
for cc=1:1:size(pa,1)
   if imag(pa(cc))==0 & pa(cc)>0
      disp(' ')
      disp(' The origin is unstable (A  has unstable real poles)')
      disp(' ')
      stop=1;
      break
   end
end
for cc=1:1:size(pa1,1)
   if imag(pa1(cc))==0 & pa1(cc)>0
      disp(' ')
      disp(' The origin is unstable (A+BC has unstable real poles)')
      disp(' ')
      stop=1;
      break
   end
end

if stop
   break
end




%%%%%%%%%%%  LMI STABILITY CONDITIONS  %%%%%%%%%%%%

tmin1=0.01;
tmin2=tmin1;

abst_init_lmi;

lmitbx_options([0 0 0 0 1]);

Q1=symmetric(n-1);

LAM1=diagonal(n-2);
lam1=LAM1*ones(n-2,1);
lam2=symmetric;

%%  Define Q0
w0=c*a*Pi;
l0=perp(w0');
F=[l0 w0'/norm(w0')];
Q0=F*[zeros(n-2) lam1; lam1' lam2]*F' + Q1;

Q0>0;
Q1>0;


dt1=(tmax1-tmin1)/(nu1-1);
timee1=tmin1:dt1:tmax1;
count=1;
for t1=timee1
	wt1=c*expm(a1*t1)*Pi;
	lt1=perp(wt1');
	Ft1=Pi'*expm(a1*t1)*Pi;
	lt1'*(Q0-Ft1'*Q1*Ft1)*lt1>0;
%	Q0-Ft1'*Q1*Ft1>0;
	count=count+1;
end

dt2=(tmax2-tmin2)/(nu2-1);
timee2=tmin2:dt2:tmax2;
count=1;
for t2=timee2
	wt2=c*expm(a*t2)*Pi;
	lt2=perp(wt2');
	Ft2=Pi'*expm(a*t2)*Pi;
	lt2'*(Q1-Ft2'*Q0*Ft2)*lt2>0;
%	Q1-Ft2'*Q0*Ft2>0;
	count=count+1;
end


lmi_mincx_tbx

Q0=value_iqc(Q0);
Q1=value_iqc(Q1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plotting min(eig(lt'*(Qi-Ft'*Q*Ft)*lt))  for tmin<t<tmax
%

ex1=10;
ex2=ex1;
clear timee1 dt1;  dt1=(tmax1+ex1)/npoi1; timee1=0:dt1:(tmax1+ex1);
count=1;
for t1=timee1
	if t1==0
	  wt1=c*a*Pi;
	  lt1=perp(wt1');
	  Ft1=eye(n-1);
	else
	  wt1=c*expm(a1*t1)*Pi;
	  lt1=perp(wt1');
	  Ft1=Pi'*expm(a1*t1)*Pi;
	end
	Rstab1(count)=min(eig(lt1'*(Q0-Ft1'*Q1*Ft1)*lt1));
%	Rstabm1(count)=max(eig(lt1'*(Q0-Ft1'*Q1*Ft1)*lt1));
%	Rstabd1(count)=min(eig(Q0-Ft1'*Q1*Ft1));
	count=count+1;
end

clear timee2 dt2;  dt2=(tmax2+ex2)/npoi2; timee2=0:dt2:(tmax2+ex2);
count=1;
for t2=timee2
	if t2==0
	  wt2=c*a*Pi;
	  lt2=perp(wt2');
	  Ft2=eye(n-1);
	else
	  wt2=c*expm(a*t2)*Pi;
	  lt2=perp(wt2');
	  Ft2=Pi'*expm(a*t2)*Pi;
	end
	Rstab2(count)=min(eig(lt2'*(Q1-Ft2'*Q0*Ft2)*lt2));
%	Rstabm2(count)=max(eig(lt2'*(Q1-Ft2'*Q0*Ft2)*lt2));
%	Rstabd2(count)=min(eig(Q1-Ft2'*Q0*Ft2));
	count=count+1;
end

%min(Rstab2)
figure(1)
clf
plot(timee1,Rstab1,'b')
hold on
%plot(timee1,Rstabm1,'b--')
%plot(timee1,Rstabd1,'b--')
plot(timee2,Rstab2,'r')
%plot(timee2,Rstabm2,'r--')
%plot(timee2,Rstabd2,'r--')
grid
title('Min(eig(R_i(t))); Want them always positive')



%min(Rstab1);
%min(Rstab2);


%%%% check if quadratic stab conditions hold for t>tmaxi

pas1=0;
pas2=0;

if min(Rstab1)>-1e-12 & max(real(pa1))<0
  pas1=check_origin(Pi'*a1*Pi,Q1,Q0,tmax1);
end

if min(Rstab2)>-1e-12 & max(real(pa))<0
  pas2=check_origin(Pi'*a*Pi,Q0,Q1,tmax2);
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
        disp(' Increase tmax1 and tmax2 and run "ofs" again')
        disp(' ')
        disp(' ')
end


