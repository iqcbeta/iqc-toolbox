function [g]=iqc_d_slope_odd_test(G)
% function [g]=iqc_d_slope_odd_test(G)
% iqc_d_slope_odd_test
%
% Parameter G can be changed to study L2 gain bounds.
% System is guaranteed stable for G<101
% Zames-Falb IQC can prove stability for G<1
% diag_slope_odd IQC can prove stability for G<101
%
% Example  proposed by A. Megretski
% Written by F. D'Amato    last modified  June 17, 1998

%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu



if nargin<1;
   G=80;  % I've changed this: ameg@mit.edu
end

s=tf([1 0],1);
e=1e-2;
M1=1/(s+1);
M2=1/(s+1+e);

fprintf(1,'System is guaranteed stable for G<101\n');
fprintf(1,'Present value of G : %f\n',G);

abst_init_iqc
lmitbx_options([0 600 -1 25 0]);
w1=signal;
w2=signal;
f=signal(2);
z1=M1*w1;
z2=M2*w2;
y=G*(z1-z2);
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

a=cell(2,2);
a{1,1}=[];a{1,2}=[100 .1];
a{2,1}=[100 .1];a{2,2}=[];

% a{1,2}=Inf;

abst_init_iqc
lmitbx_options([0 600 -1 25 0]);
w=signal(2);
f=signal(2);
z1=M1*w(1);
z2=M2*w(2);
z=[z1;z2];
y=G*(z1-z2);
v1=f(1)+y;
v2=f(2)+y;
v=[v1;v2];
[w1,xa,xb,xc,xd]=iqc_d_slope_odd(v,a);
w==w1;
g2=iqc_gain_tbx(f,z);
if isempty(g2);
   fprintf(1,'diag_slope_odd IQC: L2 gain f-->y <= infinity\n');
   g2=Inf;
else 
   fprintf(1,'diag_slope_odd IQC: L2 gain f-->y <= %f\n',g2);
end
g=[g1;g2];

% iqc_value;
% value_iqc(xa)
% value_iqc(xb)
% value_iqc(xc)
% value_iqc(xd)
