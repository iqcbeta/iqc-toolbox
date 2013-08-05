function [n,m]=size(a,b)
% function [n,m]=size(a,b)
%
% size of an "abst" object
%
% Written by ameg@mit.edu,  last modified October 13, 1997
global ABST

na=double(a);
if isempty(na)
   n=0;
   m=0;
else
   n=ABST.log(na,2);
   m=ABST.log(na,3);
end

if (nargout<2)&(nargin==1),
   n=[n m];
elseif (nargout<2)&(b==2),
   n=m;
end
