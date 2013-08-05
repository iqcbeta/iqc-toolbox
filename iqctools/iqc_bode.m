function iqc_bode(a,b)
% function iqc_bode(a,b)
% 
% plots the eigenvalues of  R(s)  vs. frequency a<s<b,
% where R(s) is the frequency weighting matrix
% of the quadratic form
%   sigma0(w,f)-sigma(Gw+Mf,w)
% sigma0(w,f)>0 is the performance specification,
% sigma(v,w)>0 is the optimal analysis IQC
%
% All eigenvalues should be non-negative
% (otherwise something's wrong), and
% the points s at which the plot is close to zero
% are expected to be the frequencies at which
% the additional IQC, if needed, should have 
% the fastest variation of their weights
% 
% the default range is set automatically
%
% Written by ameg@mit.edu

% safety checks
global ABST
A=ABST;
if ~isfield(A,'name'), 
   error('"abst" environment not defined')
end
if ~strcmp(A.name,'iqc'),
   error('This is not an "iqc" environment')
end
if ~isfield(A,'xopt'),
   error('You must run an optimizer (e.g. "iqc_gain_tbx") first')
end

E=A.E;

cst=1;
var=2;
lin=3;
lmi=4;
inp=5;
sgn=6;
vsg=7;
csg=8;
vcs=9; 
qfm=10;
iqc=11;
lnk=12;

close(gcf);

kvar=find(E.T==var);
for k=1:length(kvar),        % values of variables
   Y{k}=decs2mats(A.xopt,E.L{kvar(k)});
end

S=zeros(E.nstates+E.ninputs); % the matrix of Sigma
for j=1:E.nsimple,
   nh=0;           % C-position counter for the right (C) terms 
   ng=0;           % C-position counter for the left (B) terms
   for i=1:size(E.X{E.simples(j)},2),
      nx=E.X{E.simples(j)}(1,i); % variable number
      if nx>0,                   % variable value
         x=Y{nx}; 
      else
         x=Y{-nx}';
      end
      L=E.B{E.simples(j)}(ng+1:ng+E.X{E.simples(j)}(2,i),:)';
      R=E.gqfm(j)*E.C{E.simples(j)}(nh+1:nh+E.X{E.simples(j)}(3,i),:);
      S=S+L*x*R; 
      nh=nh+E.X{E.simples(j)}(3,i);
      ng=ng+E.X{E.simples(j)}(2,i);
   end
end
S=S+S';

c_z=E.C{E.no};          % induced L2 gain terms
c_f=E.C{E.ni};
S=A.xopt(length(A.xopt))*c_f'*c_f-c_z'*c_z-S;

sss=ss(E.ab(:,1:E.nstates), ...
   E.ab(:,E.nstates+1:E.nstates+E.ninputs), ...
   [eye(E.nstates);zeros(E.ninputs,E.nstates)], ...
   [zeros(E.nstates,E.ninputs);eye(E.ninputs)]);
if nargin==0,
   w=linspace(0,0.5/min(abs(real(eig(E.ab(:,1:E.nstates))))));
elseif nargin==1,
   w=linspace(0,a,200);
else
   w=linspace(a,b,200);
end
H=freqresp(sss,w);

nw=length(w);
y=zeros(E.ninputs,nw);
for k=1:nw,
   y(:,k)=sort(real(eig(H(:,:,k)'*S*H(:,:,k))));
end
plot(w,y);
axis([w(1) w(nw) min(y(1,:)) max(y(1,:))]);
grid





