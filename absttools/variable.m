function x=variable(N)
% function x=variable(N)
%
% creates an abstract matrix variable x with
% scalar variables placed as given by matrix N
% (interpretation same as in the LMI Control Toolbox
% except that the variable numbers are "relative",
% i.e. they will correspond to new decision variables)
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

if nargin<1,
    N=1;
end



y=abst(N);                        % convert N to "abst"

l=zeros(1,ABST.mlog);             % prepare new log entry
l(1:8)=[2 size(N) 0 5 double(y) 0 0];
z=abst_alloc(l);
x=abst(z,0);                      % let x point to this entry
