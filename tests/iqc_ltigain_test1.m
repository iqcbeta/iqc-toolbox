function iqc_ltigain_test1
% function iqc_ltigain_test1
%
% test function for "iqc_ltigain_test1", estimates the induced gain in
% v=Gw+Hf, w=d*v, d=const, -1<=d<=1
% where G=[0 g;0 0], H=[h;1]
% g(s)=(g1*s+g2)/(s^2+g3*s+g4),  h(s)=(h1*s+h2)/(s^2+h3*s+h4)
% gk,hk generated randomly
%
% result is compared to the analytically  calculated gain given by
% gamma=max{||h-g||^2+||h+g||^2} 
% where ||*|| is the H-infinity norm, gamma is gain^2.
%
% Written by ameg@mit.edu,  last modified October 13, 1997

g=rand(1,4);                   % coefficients of g(s)
h=rand(1,4);                   % coefficients of h(s)
sg=tf(g(1:2),[1 g(3:4)]);      % CS Toolbox representation of g
sh=tf(h(1:2),[1 h(3:4)]);      % CS Toolbox representation of h
gamma=1+max(norm(sg+sh,inf),norm(sg-sh,inf))^2;  % analytical gamma
G=[0 sg;0 0];
H=[sh;1];

abst_init_iqc

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sdpt3')

w=signal(2);
f=signal;
v=G*w+H*f;
w==iqc_ltigain(v); %#ok<*EQEFF>
iqc_gain_tbx(f,v)
disp(sqrt(gamma))
iqc_bode
