function x=OO(n,m)
% function x=OO(n,m)
%
% produces the abstract zero constant of size nxm
% default n=1; m=n
%
% the corresponding log entry is [1 n m 0 6 0 0 0]
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/18

if nargin==0, n=1;m=1; end
if nargin==1, m=n; end
if nargout~=1
    disp_str(3)
end

global ABST
if ~isfield(ABST,'log'),
    disp_str(12)
end

l=zeros(1,ABST.mlog);             % prepare new log entry
l(1:8)=[1 n m 0 6 0 0 0];
z=abst_alloc(l);
x=abst(z,0);                      % let x point to this entry
