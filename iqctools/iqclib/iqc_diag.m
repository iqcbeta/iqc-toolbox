function [w,X]=iqc_diag(v)
% function [w,X]=iqc_diag(v)
%
% adds iqc-terms describing the block
%
%          [d1          ]
%          |  d2        |
%  w(t) == |    d3      | * v(t)
%          |      d4    |
%          [        ... ]
%
%  where di(t) is time varying scalar with |di|<=1
%
% Created cykao@mit.edu, last modified on Feb 16, 1998

m=size(v,1);  
w=signal(m);
X=diagonal(m);
X>0;
v'*X*v-w'*X*w>0;
