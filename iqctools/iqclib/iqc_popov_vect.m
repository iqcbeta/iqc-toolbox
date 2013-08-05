function [w,Lambda]=iqc_popov_vect(v,sign)
% function [w,Lambda]=iqc_popov_vect(v,sign)
%
% Popov IQC for sector nonlinearities and parametric
% uncertainty (any size of v).
%
% w'*Lambda*(dv/dt)>0, 
%
% where w=phi(v) in the case of nonlinearity and w=delta*v 
% for the case with parametric uncertainty
%
% Here Lambda is symmetric nxn matrix with
%
% Lambda>0 if sign='+'
% Lambda<0 if sign='-'
% Lambda unconstrained if sign='0', which is default
%

%
% Work by Ulf Jonsson Nov 1997
%

if nargin<1
   error('input must be specified')
end
if nargin<2
  sign='0';
end
n=size(v,1);
u=derivative(v);
w=signal(n);
Lambda=symmetric(n);
if sign=='+'
  Lambda>0;
elseif sign=='-'
  Lambda<0;
end
w'*Lambda*u>0;


