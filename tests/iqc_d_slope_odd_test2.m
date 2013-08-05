% iqc_d_slope_odd_test2
% function [g]=iqc_d_slope_odd_test2(P,C,option5);%
%
% Written by fdmato@lids.mit.edu    last modified  June 22, 1998

function [g]=iqc_d_slope_odd_test2(P,C,option5);

s=tf([1 0],1);
if nargin<3;option5=1;end
if nargin<2|isempty(C);
  C=-1;
end
if nargin<1|isempty(P);
  P=1;
end

d=1e-2;
M1=1/(s+1);
M2=1/(s+1+d);

abst_init_iqc
lmitbx_options([0 600 -1 25 option5]);
  w1=signal;
  w2=signal;
  f=signal;
  y=signal;

  v1=M1*(f+C*y);
  v2=M2*(f+C*y);
  w1==iqc_monotonic(v1,1,[1 100]);
  w2==iqc_monotonic(v2,1,[1 100]);
  y==P*(w1-w2);
  g1=iqc_gain_tbx(f,[v1;v2]);
if isempty(g1);
 fprintf(1,'Zames-Falb IQC: L2 gain f-->y <= infinity\n');
 g1=Inf;
else 
 fprintf(1,'Zames-Falb IQC: L2 gain f-->y <= %g\n',g1);
end 

a={100,Inf;Inf,100};

abst_init_iqc
lmitbx_options([0 600 -1 25 option5]);	
  w=signal(2);
  f=signal;
  y=signal;

  v1=M1*(f+C*y);
  v2=M2*(f+C*y);
  [waux,xa,xb,xc,xd,dd]=iqc_d_slope_odd([v1;v2],a);
  w==waux;
  y==P*(w(1)-w(2));
  g2=iqc_gain_tbx(f,y);
if isempty(g2);
 fprintf(1,'d_slope_odd IQC: L2 gain f-->y <= infinity\n');
 g2=Inf;
else 
 fprintf(1,'d_slope_odd IQC: L2 gain f-->y <= %f\n',g2);
end 
g=[g1;g2];


