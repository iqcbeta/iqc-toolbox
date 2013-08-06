function [w,X,Y,Z,Lambda]=iqc_polytope_stvp(v,Delta,Omega,Lambdastruc)
%function [w,X,Y,Z,Lambda]=iqc_polytope_stvp(v,Delta,Omega,Lambdastruc);
%
% Defines the IQC for a slowly time-varying polytopic uncertainty,
% i.e., w(t)=Delta(t)v(t), where Delta(t) takes values in the
% polytope C1={Delta1,...,DeltaN} and d/dt(Delta) takes values in
% the polytope C2={Omega1,...,OmegaN}. We assume that 0 is in C1
% and C2. The IQC is
%
%  [v]'[Z   Y][v]
%  [ ]*[     ][ ]+2*w'*Lambda*(dv/dt)>0
%  [w] [Y' -X][w]
%
%  where X=X'>0, Z=Z', Lambda has structure Lambdastruc, and
%
%  [I ]'[Z   Y][I ]
%  [  ]*[     ][  ]-1/2(Lambda'*Oj+Oj'*Lambda)>0, for all i,j
%  [Di] [Y' -X][Di]
%
%  Di'*Lambda=Lambda'*Di,  i=1,...,N  (*)
%
%  where Di=Delta{i}, Oj=Omega{j}
%
%  Inputs: Delta and Omega are cellarays containing the vertices
%          of the polytopes C1 and C2.
%          Lambdastruc is a structure matrix such that (*) holds.
%

%
% Work by Ulf Jonsson, Nov 1997
%

if nargin<4,
    disp_str(4,'Four')
end
N1=size(Delta,2);
N2=size(Omega,2);
nx=size(Delta{1},1);
nz=size(Delta{1},2);
if size(Omega{1},1) ~= nx
    error('Different size polytopes C1 and C2')
end
if size(Omega{1},2) ~= nz
    error('Different size polytopes C1 and C2')
end
if size(v,1)~=nz
    error('Wrong size of signal v')
end

global ABST
switch ABST.systemtype
    case 'continuous'
        w=signal(nx);
        u=derivative(v);
        X=symmetric(nx);
        X>0; %#ok<*VUNUS>
        Z=symmetric(nz);
        Y=rectangular(nz,nx);
        Lambda=variable(Lambdastruc);
        for i=1:N1
            Di=Delta{i};
            for j=1:N2
                Oj=Omega{j};
                Z+Y*Di+Di'*Y'-Di'*X*Di-(1/2)*(Lambda'*Oj+Oj'*Lambda)>0;
            end
        end
        v'*Z*v+v'*Y*w+w'*Y'*v-w'*X*w+(2*w)'*Lambda*u>0;
    case 'discrete'
        disp_str(70,'iqc_polytope_stvp')
end
