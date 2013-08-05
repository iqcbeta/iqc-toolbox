function g=derivative(f)
% function g=derivative(f)
%
% g=df/dt
%
% Written by ameg@mit.edu,  last modified October 13, 1997

global ABST

if ABST.log(double(f),1)~=6 & ABST.log(double(f),1)~=5
   error('only a signal or an input could be differentiated')
end

l=zeros(1,ABST.mlog);             % prepare new log entry
l(1:8)=[6 size(f,1) 1 -3 0 double(f) 0 0];
z=abst_alloc(l);
g=abst(z,0);                      % let x point to this entry