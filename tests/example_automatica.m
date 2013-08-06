
%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu


% --- data ---(from Kothare's PHD thesis)
s=tf([1 0],1);
G=10/(100*s+1)*[4 -5 ; -3 4];
Q=(100*s+1)/10/(20*s+1)*[ 4 5;3 4];
% ff=2.5*(s+1)*eye(2,2);
Gtilde=10/(100*s+1)*[4          -5/(0.1*s+1);
		    -3/(.1*s+1)      4];
Q1=2.5*(s+1)/(20*s+1)/(.1*s+1)*[1.6*s+1   2*s;
		                1.2*s    1.6*s+1];
Q2=1/(100*s+1)*[99               -125*(s+1)/(.1*s+1);
               -75*(s+1)/(.1*s+1)    99];
	       
% --- initialize IQC environment ---
abst_init_iqc;

w=signal(2);
f=signal(2);

z = G*(w+f);
e = z-Gtilde*w;
v = -Q1*e-Q2*w;


% --- Zames-Falb approach ---                       --- gamma - ---
% location={[Inf 0.3] [];[] [Inf 0.3]};              %      4.1323 
%location={[Inf best] [];[] [Inf best]};             %      4.1323 

% --- quadratic lyapunov function ---
location={Inf Inf;Inf Inf};                      
[waux,xa,xb,xc,xd,dd]=iqc_d_slope_odd(v,location,0,1);   %    7.8156

                                                % LMI gave  7.8523
                                                % (close, but does not
                                                %  verify frequency
                                                % domain constraints)
% --- new results (generalized Zames-Falb IQC) ---
location={[Inf .3] [ Inf .3];[Inf .3] [Inf .3]};    %      3.2327

[waux,xa,xb,xc,xd,dd]=iqc_d_slope_odd(v,location,0,1);

w==waux;

g=iqc_gain_tbx(f,w);


iqc_value
xxa=value_iqc(xa);
xxb=value_iqc(xb);
if ~isa(xc,'abst');
  xxc=[];
else
  xxc=value_iqc(xc);
end
if ~isa(xd,'abst');
  xxd=[];
else
  xxd=value_iqc(xd);
end

ddd=value_iqc(dd);

% recover multipliers G and H 
[GG,HH]=getGH_d_slope_odd(xxa,xxb,xxc,xxd,ddd,location);
