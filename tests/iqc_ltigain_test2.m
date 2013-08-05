function iqc_ltigain_test2(cs)
% function iqc_ltigain_test2(cs)
%
% estimates the gain from v to y in
% y=gw+hv, w=d*v, d=const, -1<=d<=1
% where
% g(s)=(g1*s+g2)/(s^2+g3*s+g4),  h(s)=(h1*s+h2)/(s^2+h3*s+h4)
% gk,hk generated randomly
%
% result is compared to the analytically  calculated gain given by
% gamma=1+max{||h-g||^2,||h+g||^2}
% where ||*|| is the H-infinity norm, gamma is gain^2.
%
% Written by ameg@mit.edu,  last modified 5.18.97
% Modified by rantzer 8.28.97
% modified by ameg October 13, 1997
if nargin>0,
   rand('seed',cs);
end

cg=rand(1,4);                   % coefficients of g(s)
ch=rand(1,4);                   % coefficients of h(s)
g=tf(cg(1:2),[1 cg(3:4)]);      % CS Toolbox representation of g
h=tf(ch(1:2),[1 ch(3:4)]);      % CS Toolbox representation of h
gamma=max(norm(g+h,inf),norm(g-h,inf));  % analytical gamma

abst_init_iqc
w=signal;
v=signal;
y=g*w+h*v;
w==iqc_ltigain(v);
iqc_gain_tbx(v,y)
disp(gamma)
iqc_bode