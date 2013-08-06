function x=skew(n)
% function x=skew(n)
%
% defines an "abst" object which is a
% skew symmetric matrix variable
% of size n by n
%
% the corresponding log entry will be
% [2 n n 0 4 0 0 0]
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/18

if nargin==0, n=2; end
if n<2,
    disp_str(7,'skew','one')
end
if nargout~=1
    disp_str(3)
end

global ABST
if ~isfield(ABST,'log'),
    disp_str(12)
end

l=zeros(1,ABST.mlog);             % prepare new log entry
l(1:8)=[2 n n 0 4 0 0 0];
z=abst_alloc(l);
x=abst(z,0);                      % let x point to this entry
