function x=II(n)
% function x=II(n)
%
% produces the abstract identity constant of size nxn
% default n=1
%
% the corresponding log entry is [1 n n 0 7 0 0 0]
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/18

if nargin==0, n=1; end
if nargout~=1
    disp_str(3)
end

global ABST
if ~isfield(ABST,'log'),
    disp_str(12)
end

l=zeros(1,ABST.mlog);             % prepare new log entry
l(1:8)=[1 n n 0 7 0 0 0];
z=abst_alloc(l);
x=abst(z,0);                      % let x point to this entry
