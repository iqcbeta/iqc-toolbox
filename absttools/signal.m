function x=signal(n)
% function x=signal(n)
%
% defines "abst" input vector signal of size n
%
% the corresponding log entry will be
% [5 n 1 0 1 0 0 0]
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/18

if nargout~=1
    disp_str(3)
end

global ABST
if ~isfield(ABST,'log'),
    disp_str(12)
end

if nargin<1, n=1; end

l=zeros(1,ABST.mlog);             % prepare new log entry
l(1:8)=[5 n 1 0 8 0 0 0];
z=abst_alloc(l);
x=abst(z,0);                      % let x point to this entry
