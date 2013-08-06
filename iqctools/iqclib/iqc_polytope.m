function [w,X,Y,Z]=iqc_polytope(v,C)
%function [w,X,Y,Z]=iqc_polytope(v,C)
% Defines the IQC for a polytopic uncertainty, i.e.,
% w(t)=Delta(t)v(t), where Delta(t) takes values in the
% polytope C={Delta1,...,DeltaN}. We assume 0 is in C
%
% The IQC is
%
%  [v]'[Z   Y][v]
%  [ ]*[     ][ ]>0
%  [w] [Y' -X][w]
%
%  where X=X'>0, Z=Z', and
%
%  [    I   ]'[Z   Y][   I    ]
%  [        ]*[     ][        ]>0, for all i
%  [Delta{i}] [Y' -X][Delta{i}]
%
%  Input: C is a cellaray containing the vertices of
%         the polytope
%
% Last modified by cmj on 2013/5/4

if nargin<2
    disp_str(4,'Two')
end

N=size(C,2);
nx=size(C{1},1);
nz=size(C{1},2);
if size(v,1)~=nz,
    error('size of v is not consistent with size of polytope')
end
w=signal(nx);
X=symmetric(nx);
X>0; %#ok<*VUNUS>
Z=symmetric(nz);
Y=rectangular(nz,nx);
M=[Z Y;Y' -X];
Id=eye(nz);
for k=1:N
    [Id;C{k}]'*M*[Id;C{k}]>0;
end
[v;w]'*(M*[v;w])>0;
