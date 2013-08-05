function [Q]=hys_globsta(G,d,nu);
%HYS_GLOBSTA  Global stability of limit cycles of RFS
%
%   HYS_GLOBSTA(SYS,D,NU)  proves global stability of the limit
%   cycle of the system SYS in hysteresis relay feedback.  It
%   calculates the Poincare' map x_k+1=F(t)x_k and finds a Q>0 such
%   that the Poincare' map is globally quadratically stable, i.e.,
%   finds Q>0 such that
%   		alpha(t)-tau(t)beta(t)>0   for all tmin<t<tmax
%   where alpha(t)=Q-F'(t)QF(t) and tmin and tmax are, respectively,
%   the minimum and maximum times between any two switches of any
%   trajectory after a long time.
%
%   SYS is an LTI system SYS (created with either TF or SS).  The system 
%   must be stable and C*inv(A)*B < 0.
%
%   D is the hysteresis parameter.  It needs to satisfy 0<D<-C*inv(A)*B.
%   It can be 0 if CB<0.  By default is set to -C*inv(A)*B/2.
%
%   NU is the number of LMIs used to sample the interval tmin<t<tmax.
%   By default is set to 40.
%
%   For example
%
%	s=tf([1 0],1);
%	G=1/((s+1)*(s+1)*(s+1));
%	hys_globsta(G,0.5,40)
%
%   Version: 1.0    Last modified: Aug 9, 1999

%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu



%clear
%s=tf([1 0],1);
%G=1/((s+1)*(s+1)*(s+1));
%epsl=0.0001;
%G=tf([-epsl 1-3*epsl 12-4*epsl],[1 3 3 4]);
%	%%Make alpha(0)=0 if zeroat0=1
%tmax=1;
%	%%number of LMIs for main interval
%nu=50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  INITIALIZATION OF VARIABLES  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ini



%%%%%%%%%%  Default variable nu  %%%%%%%%%%%%%%%%
ttt=which('nu');
if isempty(ttt)
	nu=40;
end
clear ttt


% bounds on tmin and tmax

tmin=bound_tmin(a,b,c,d);
tmax=bound_tmax(a,b,c,d,tmin);

disp(['Stability condition will be sampled in the interval [' num2str(tmin) ',' num2str(tmax) ']'])
disp(' ')



%%%%%%%%%   at t == 0   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q0=c*a*Pi;
p0=-c*a*(xstar+ainv*b);       %%%  = -(c*b+c*a*xstar);
l0=q0'*p0/(q0*q0');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Main LMI's  definition  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


abst_init_lmi;

lmitbx_options([0 0 0 0 1]);

Q=symmetric(n-1);
Q>0;
tau=diagonal(nu);
tau>0;

dt=(tmax-tmin)/(nu-1);
timee=tmin:dt:tmax;
count=1;
for t=timee
	vt=(expm(a*t)-expm(a*tstar))*barxstar;
	qt=c*expm(a*t)*Pi;
	pt=-c*vt;
	lttemp=qt'*pt/(qt*qt');
	sttemp=p0*(l0/norm(l0)*norm(lttemp)-lttemp/norm(lttemp)*norm(l0));
	sttemp=sttemp/norm(sttemp);
	betat=lttemp*sttemp'+sttemp*lttemp';
	ltime(:,count)=lttemp;
	st(:,count)=sttemp;
	Ft=Pi'*(vt*c/(c*vt)-eye(n))*expm(a*t)*Pi;
	Q-Ft'*Q*Ft-tau(count,count)*betat>0;
	count=count+1;
end

%%%%%%%%%%%  imposing LMI at t*

Ftstar=Pi'*(v*c/(c*v)-eye(n))*expm(a*tstar)*Pi;
qtstar=c*expm(a*tstar)*Pi;
ztstar=perp(qtstar');

ztstar'*(Q-Ftstar'*Q*Ftstar)*ztstar>0;

lmi_mincx_tbx
%lmi_mincx_tbx(trace(Q))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  End of Main LMI's  definition  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Q=value_iqc(Q);
tau=value_iqc(tau);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plotting max(eig(Q-Ft'*Q*Ft-tau*betat)) for tmin<t<tmax
%

count=1;
for t=timee
	betat=ltime(:,count)*st(:,count)'+st(:,count)*ltime(:,count)';
	vt=(expm(a*t)-expm(a*tstar))*barxstar;
	Ft=Pi'*(vt*c/(c*vt)-eye(n))*expm(a*t)*Pi;
	Rstab(count)=min(eig(Q-Ft'*Q*Ft-tau(count,count)*betat));
%	Rstabmax(count)=max(eig(Q-Ft'*Q*Ft-tau(count,count)*betat));
	Rstabalt(count)=min(eig(Q-Ft'*Q*Ft));
	count=count+1;
end

figure(1)
plot(timee,Rstab)
grid
title('min(eig(Q-Ft''*Q*Ft-tau*betat)).  Want it to be always >0')
%plot(timee,Rstabmax,'r')
%plot(timee,Rstabalt,'k')

positive=min(Rstab);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  FINAL CONCLUSIONS  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%test if  Q-Ft'*Q*Ft-taut*betat>0 for tmin<t<tmax

% if feasible then limit cycle is globally asymp. stable.

if (positive>0)
	disp(' ')
	disp(' ')
	disp('   THE LIMIT CYCLE IS GLOBALLY ASYMPTOTICALLY STABLE')
	disp(' ')
	disp(' ')
else
	disp(' ')
	disp(' ')
	disp('   The limit cycle could not be proven to be globally stable')
	disp(' ')
	disp(' ')
        Q=[];
end



