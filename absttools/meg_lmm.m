function [y,x,nsteps,axb]=meg_lmm(a,b,ymin,x0,maxstep,tol)
% function [y,x,nsteps,axb]=meg_lmm(a,b,ymin,x0,maxstep,tol)
% 
% solves the linear program   max_j{a(:,j)'*x+b(j)}->min_x
% does not care to minimize below ymin, defauilt ymin=-1e8,
% with the initial guess x0 (column vector, default x0=0)
% maximal # of iterations maxstep, default maxstep=10000,
% tolerance parameter tol, default tol=1e-8,
% 
%
% outputs: y         =   min. value
%          x         =   optimal vector of decision variables
%          nsteps    =   # of iterations performed
%          axb       =   a(:,j)'*x+b(j)
%
% this M-file calls a mex-routine lmm written by ameg@mit.edu

% first check the inputs and set default values
if nargin<2, 
   error('Two first arguments must be  defined'); 
end
n=size(a,1);  % # of decision variables
m=size(a,2);  % # of linear terms
if ~isequal([1 m],size(b)), 
   b=b';
   if ~isequal([1 m],size(b)),
      error('Improper second argument size');
   end
end
if nargin<3,
   ymin=-1e8;
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

% now apply lmm
[y,x,nsteps,axb]=lmm(a,b,x0,ymin,tol,maxstep);

