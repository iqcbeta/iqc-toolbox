function y=abst_trace(x)
% function y=abst_trace(x)
%
% Written by ameg@mit.edu,  last modified October 13, 1997
% Last modified by cmj on 2013/4/18

n=size(x,1);
if n~=size(x,2),
    disp_str(34)
end

y=x(1,1);
for i=2:n,
    y=y+x(i,i);
end
