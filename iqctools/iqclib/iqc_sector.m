function [w,x]=iqc_sector(v,alpha,beta)
% function [w,x]=iqc_sector(v,alpha,beta)
%
% Defines the IQC for a sector bounded scalar nonlinearity, i.e.
% 
%       x(w-alpha*v)*(beta*v-w)>0,
%
% where x>0 and w=phi(v).
%
% Default values: alpha=0
%                 beta=1
%
% Work by ulfj@mit.edu Oct 28, 1997.
% Last modified by cmj on 2013/5/4

if nargin<1
   disp_str(1)
end
if size(v,1)~=1
   disp_str(8,'scalar')
end
if nargin<2
    alpha=0; 
end
if nargin<3
    beta=1;
end

x=symmetric;
x>0;
w=signal;
(w-alpha*v)'*x*(beta*v-w)>0;
