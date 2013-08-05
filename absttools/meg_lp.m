function [y,x,nsteps]=meg_lp(A,B,C,x0,maxstep,tol)
% function [y,x,nsteps]=meg_lp(A,B,C,x0,maxstep,tol)
% 
% solves the linear program   Ax>B, y=Cx->min_x
% with the initial guess x0 (column vector, default x0=0)
% maximal # of iterations maxstep, default maxstep=10000,
% tolerance parameter tol, default tol=1e-8,
% 
%
% outputs: y         =   min. value
%          x         =   optimal vector of decision variables
%          nsteps    =   # of iterations performed
%
% this M-file calls a mex-routine lmm written by ameg@mit.edu

% first check the inputs and set default values
if nargin<3, 
   error('Three first arguments must be  defined'); 
end
n=size(A,2);  % # of decision variables
m=size(A,1);  % # of linear constraints
if ~isequal([m 1],size(B)), 
   error('Improper second argument size');
end
if ~isequal([1 n],size(C)),
   error('Improper third argument size');
end
if nargin<4,
   x0=zeros(n,1);
elseif ~isequal([n 1],size(x0)),
   error('Improper third argument size');
end
if nargin<5,
   maxstep=10000;
end
if nargin<6,
   tol=0.00000001;
end

% finding a feasible solution
[y1,x1,nsteps1,axb]=lmm(-A',B',x0,-1,tol,maxstep);
if ~(y1<0),
   y=[];
   nsteps=nsteps1;
   x=x1;
   return;
end

% preparing for the second lmm run
a=[-A' C'];
b=[zeros(1,m) 1];
for k=1:m,
   a(:,k)=a(:,k)*(-1/axb(k));
end

[y2,x2,nsteps2,axb]=lmm(a,b,zeros(n,1),-1,tol,maxstep);
if y2<=0,
   y=-Inf;
   x=x2;
else
   y=1-1/y2+C*x1;
   x=x1+x2/y2;
end
nsteps=nsteps1+nsteps2;



