% iqc_d_slope_odd_test1
%
% Parameter C can be changed to study L2 gain bounds.
% System is guaranteed stable for |C| < 101
% Zames-Falb IQC can prove stability for |C|<1
% diag_slope_odd IQC can prove stability for -793 < C < 101
%
% Example proposed by ameg@mit.edu
% Written fdamato@lids.mit.edu       last modified  June 17, 1998

function [g]=iqc_d_slope_odd_test(C)

if nargin<1;
    C=100;
end

s=tf([1 0],1);
e=1e-2;
M1=1/(s+1);
M2=1/(s+1+e);

fprintf(1,'SGT guarantees stability for |C|<101\n');
fprintf(1,'Present value of C : %f\n',C);

abst_init_iqc
lmitbx_options([0 600 -1 25 1]);

w1=signal;
w2=signal;
f=signal(2);

z1=M1*w1;
z2=M2*w2;
y=C*(z1-z2);
v1=f(1)+y;
v2=f(2)+y;
z=[z1;z2];
w1==iqc_monotonic(v1,[0.1 3],1);
w2==iqc_monotonic(v2,[0.1 3],1);
g1=iqc_gain_tbx(f,z);
if isempty(g1);
    fprintf(1,'Zames-Falb IQCs: L2 gain f-->y <= infinity\n');
    g1=Inf;
else
    fprintf(1,'Zames-Falb IQCs: L2 gain f-->y <= %g\n',g1);
end

a={[] Inf; Inf []};

abst_init_iqc
% lmitbx_options([0 600 -1 25 1]);

w=signal(2);
f=signal(2);

z1=M1*w(1);
z2=M2*w(2);
z=[z1;z2];
y=C*(z1-z2);
v1=f(1)+y;
v2=f(2)+y;
v=[v1;v2];
[waux,xa,xb,xc,xd,dd]=iqc_d_slope_odd(v,a);
w==waux;
g2=iqc_gain_tbx(f,z);

if isempty(g2);
    fprintf(1,'d_slope_odd IQC: L2 gain f-->y <= infinity\n');
    g2=Inf;
else
    fprintf(1,'d_slope_odd IQC: L2 gain f-->y <= %f\n',g2);
end
g=[g1;g2];





