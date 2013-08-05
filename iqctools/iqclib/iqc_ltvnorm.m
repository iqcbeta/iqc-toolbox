function [w,x]=iqc_ltvnorm(v,n,gain)
% function [w,x]=iqc_ltvnorm(v,n,gain)
%
% defines w (size n) as related to v by the L2 norm bound
% ||w||<gain*||v||
% 
% x>0 is the salar multiplier of the IQC gain^2*v'*x*v>w'*x*w
%
% Written by ameg@mit.edu,  last modified October 13, 1997
if nargin<1,
   error('input must be specified')
end
if nargin<2, n=size(v,1); end
if nargin<3, gain=1; end

x=symmetric;
w=signal(n);
x>0;
if gain==1,
   v'*x*v>w'*x*w;
else
   v'*(gain^2)*x*v>w'*x*w;
end
