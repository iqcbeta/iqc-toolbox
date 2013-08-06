function x=rectangular(n,m)
% function x=rectangular(n,m)
%
% defines an "abst" object which is a
% full rectangular matrix variable
% of size n by m
%
% the corresponding log entry will be
% [2 n m 0 2 0 0 0]
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/18

if nargin<1, n=1; m=1; end
if nargin<2, m=n; end
if nargout~=1
    disp_str(3)
end

global ABST
if ~isfield(ABST,'log'),
    disp_str(12)
end

l=zeros(1,ABST.mlog);             % prepare new log entry
l(1:8)=[2 n m 0 2 0 0 0];
z=abst_alloc(l);
x=abst(z,0);                      % let x point to this entry
