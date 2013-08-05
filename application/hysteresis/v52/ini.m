%
% System Definition, existence, and local stability
%


%clear

figure(1)
clf
figure(2)
clf
hold on
grid


[a,b,c]=ssdata(G);


%%% Testing the stability of the system %%%

if max(real(eig(a)))>=0
	disp(' ')
	disp('THE OPEN LOOP SYSTEM IS NOT STABLE')
	disp(' ')
	disp(' ')
	break
end	

ainv=inv(a);


%%% Testing if the equilibrium point is on the right side %%%

if c*ainv*b>=0
	disp(' ')
	disp('	THE SYSTEM HAS EQUILIBRIUM POINTS')
%	disp('	THE SYSTEM CANNOT HAVE GLOBALLY STABLE LIMIT CYCLES')
	disp('		NEED TO HAVE  C*inv(A)*B < 0')
	disp(' ')
	disp(' ')
	break
else
	disp(' ')
	disp('- System has no equilibrium points')	
end


%%%%%%%%%%  Default variable d  %%%%%%%%%%%%%%%%
ttt=which('d');
if isempty(ttt)
        d=-c*ainv*b*0.05;
end
clear ttt


%%% Checking if d is good, i.e., if d>0 or d>=0 and d<-C*inv(A)*B %%%


if c*b<0
   if ((d>= -c*ainv*b) | (d<0))
	disp(' ')
	disp(['	BAD CHOICE OF d.  d must be in [0,' num2str(- c*ainv*b) ')'])    %)
	disp(' ')
	disp(' ')
	break
   else
	disp(['- Correct choice of d between [0,' num2str(- c*ainv*b) ')'])
   end
else	  % c*b>=0
   if ((d>= -c*ainv*b) | (d<=0))
	disp(' ')
	disp(['	BAD CHOICE OF d.  d must be in (0,' num2str(- c*ainv*b) ')'])    %)
	disp(' ')
	disp(' ')
	break
   else
	disp(['- Correct choice of d between (0,' num2str(- c*ainv*b) ')'])
   end
end


% bounds on tmin and tmax

tmin=bound_tmin(a,b,c,d);
tmax=bound_tmax(a,b,c,d,tmin);


n=size(a,1);

%%% String variables %%%
A=mat2str(a);
B=mat2str(b);
C=mat2str(c);
Ainv=mat2str(ainv);
N=mat2str(n);
D=mat2str(d);

%%%%  Finding limit cycle(s)  %%%%%%%%%%%%%%%%%%

g=inline([C '*inv((expm(' A '*t)+eye(' N ')))*(expm(' A '*t)-eye(' N '))*' Ainv '*' B '-' D]);

%
% Finding zeros of g
%

time_zero_g=ezero(g,0,tmax,0.05);


%%%%  TESTING THE LIMIT CYCLES to see if they are >-d on [0,tstar] %%%%

num_lc=0;
time_tstar=[];
time_not_tstar=[];
for count=1:size(time_zero_g,2)
        tstar_c=time_zero_g(count);
        xstar_c=inv((expm(a*tstar_c)+eye(n)))*(expm(a*tstar_c)-eye(n))*ainv*b;
        ok_lc=1;
        for tii=0.01:0.1:(tstar_c-0.01)
	        if (c*expm(a*tii)*(xstar_c-ainv*b)+c*ainv*b+d)<0
			ok_lc=0;
        	        break
		end
	end
	if ok_lc
		num_lc=num_lc+1;
                time_tstar=[time_tstar tstar_c];
        else
                time_not_tstar=[time_not_tstar tstar_c];
	end
end

%time_zero_g
%time_tstar
%time_not_tstar


%%%%%%%   At this point I know how many limit cycles exist and where they are.
%%%%%%%   Let's plot their trajectories and also the g function.

figure(2)
title('Checking existence of symmetric unimodal limit cycles')

if num_lc==0 
        fplot(g,[0;tmax])
        disp(' ')
        disp(' ')
        disp('THE SYSTEM HAS NO SYMMETRIC UNIMODAL LIMIT CYCLES')
        disp(' ')
        disp(' ')
        break
end



%%%%  if the program is here is because there exists zeros of g


%%% plotting y+d of time_not_tstar

for count=1:size(time_not_tstar,2)
     tstar=time_not_tstar(count);
     xstar=inv((expm(a*tstar)+eye(n)))*(expm(a*tstar)-eye(n))*ainv*b;
     XSTAR=mat2str(xstar);
     y=inline([C '*expm(' A ' *t)*(' XSTAR '-' Ainv '*' B ')+' C '*' Ainv '*' B '+' D]);
     fplot(y,[0,tstar],'--g')
end

%%% plotting y+d of time_tstar

for count=1:size(time_tstar,2)
     tstar=time_tstar(count);
     xstar=inv((expm(a*tstar)+eye(n)))*(expm(a*tstar)-eye(n))*ainv*b;
     XSTAR=mat2str(xstar);
     y=inline([C '*expm(' A ' *t)*(' XSTAR '-' Ainv '*' B ')+' C '*' Ainv '*' B '+' D]);
     fplot(y,[0,tstar],'k')
end


%%%  ploting the g function
time_scale=max([time_tstar time_not_tstar]);
fplot(g,[0;1.5*time_scale])


%legend('y for x(0)=x* (no switch)','relay hysteresis','relay hysteresis','g function')


if num_lc>1
	disp(' ')
	disp('THERE EXISTS MORE THAN ONE LIMIT CYCLE')
%        disp(' Limit cycles half period:')
%        time_tstar
	disp(' ')
	break
end

disp(['- The system has one symmetric unimodal limit cycle with period 2* ' num2str(tstar)])

barxstar=xstar - ainv*b;

Pi= perp(c');


%%% Checking local stability of the Poincare' map %%%

v=-a*xstar-b;
w=(eye(n)-v*c/(c*v))*expm(a*tstar);
eg=eig(w);

m=max(abs(eg));

if m>=1
	disp('The limit cycle is locally unstable')
	disp(' ')
	break
else
	disp('- Limit cycle is locally stable')
	disp(' ')
end


%%% ckech c*vt


BARXSTAR=mat2str(barxstar);
TSTAR=mat2str(tstar);      
cvt=inline([C '*(expm(' A '*t)-expm(' A '*' TSTAR '))*' BARXSTAR ]);

time_zeros_cvt=ezero(cvt,tmin,tmax,0.05);

if size(time_zeros_cvt,2)>1
    disp(['  **Warning**   c*vt = 0 at   ' num2str(time_zeros_cvt)])
    disp(' ')
end


%figure(3)
%fplot(cvt,[0,10])
%grid


clear g dt t e antes gfun time_zero_aprox count time_zero_g num_lc time_tstar tstar_c xstar_c ok_lc tii XSTAR BARXSTAR TSTAR Ainv D N A B C time_zeros_cvt y w eg m time_scale time_not_tstar cvt

