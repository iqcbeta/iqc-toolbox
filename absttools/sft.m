function y=sft(x,n)
% function y=sft(x,n)
%
% fits  string x to length n, 
% padding with empty spaces on the right, if necessary
% 
% default n=10
%
% Written by ameg@mit.edu,  last modified October 13, 1997
if nargin<2, 
 n=10; 
end

if ischar(x),
 z=x;
else
 z=num2str(x);
end
 

if length(z)>=n,
 y=z(1:n);
else
 y=[z blanks(n-length(z))];
end

 
