function y=abst_trace(x)
% function y=abst_trace(x)
%
% Written by ameg@mit.edu,  last modified October 13, 1997
n=size(x,1);
if n~=size(x,2),
   error('argument of "trace" must be square')
end

y=x(1,1);
for i=2:n,    
    y=y+x(i,i);
end
