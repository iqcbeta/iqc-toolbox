function x=symmetric(n)
% function x=symmetric(n)
%
% defines an "abst" object which is a 
% full symmetric matrix variable
% of size n by n
%
% the corresponding log entry will be
% [2 n n 0 1 0 0 0]
%
% Written by ameg@mit.edu,  last modified October 13, 1997
global ABST
if ~isfield(ABST,'log'),
   error('"abst" environment not initialized')
end

if nargin<1, n=1; end

l=zeros(1,ABST.mlog);             % prepare new log entry
l(1:8)=[2 n n 0 1 0 0 0];
z=abst_alloc(l);
x=abst(z,0);                      % let x point to this entry
