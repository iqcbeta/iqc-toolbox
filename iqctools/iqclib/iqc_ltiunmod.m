function [w,x]=iqc_ltiunmod(v,a,n,gain)
% Use: [w,x]=iqc_ltiunmod(v,a,n,gain)
% 
% Defines w (size n) which relates to v by 
%
%  w= Delta(s) v     
%  
% where Delta is an (uncertain) lti system 
% with || Delta(s)||< gain.
% 
% Arguments 
%          a: (-) poles of multiplier
%          n:  output size
%          gain : norm bound
%
% Written by F. D'Amato,   last modified: June 11, 1998

ni=size(v,1);
if nargin < 4; gain=1; end
if nargin < 3; n=ni; end
if nargin < 2; a=1; end

k2=gain^2;
w=signal(n);
na=length(a);
x=rectangular(na,1);
x0=rectangular; 

s=tf([1 0],1);

p=x0;
v'*k2*x0*v > w'*x0*w;

for ndx=1:na;
   G=a(ndx)/(s+a(ndx));   
   vg=G*v;
   wg=G*w;
   p=p+x(ndx)*G;
   v'*k2*(x(ndx)*vg) > w'*(x(ndx)*wg);
end

p>0;
