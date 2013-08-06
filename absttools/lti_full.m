function [y,x]=lti_full(num,den,u)
% function [y,x]=lti_full(num,den,u)
%
% this function helps get access to states of lti transformations
% of abst signals (inp and sgn)
% here y=[num(s)/den(s)]*u,
% x is the cell array with n+1 elements
% x{k}=(d/dt)^{k-1}*[1/den(s)]*u.
% this works for VECTOR u
%
% Written by ameg@mit.edu,  last modified Dec. 07, 1997

n=length(den)-1;
num=[zeros(1,n+1-length(num)) num];
v=tf(1,den)*u;
x{1}=v;
y=num(n+1)*v;
for k=1:n,
    v=derivative(v);
    y=y+num(n+1-k)*v;
    x{k+1}=v; %#ok<*AGROW>
end
