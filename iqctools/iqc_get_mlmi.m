function [P,A,B,Sigma,MainLMI]=iqc_get_mlmi
% function [P,A,B,Sigma,MainLMI]=iqc_get_mlmi
% 
% extract P, A, B, and Sigma matrices of the
% main LMI of the optimization problem, i.e.
%
% MainLMI == [PA+A'P, PB; B'P, 0]+Sigma < 0
%
% written by cykao@mit.edu,  last modified Oct. 21 1998

% safety checks
global ABST
A=ABST;
if ~isfield(A,'name'),
   error('"abst" environment not defined')
end
if ~strcmp(A.name,'iqc'),
   error('This is not an "iqc" environment')
end
if ~isfield(A,'P'),
   error('You must run optimizer "iqc_gain_tbx" first')
end

% extracting P, A, B matrices
P=ABST.P;
AB=ABST.E.ab;
n=size(AB,1);
m=size(AB,2);
A=AB(:,1:n);
B=AB(:,n+1:m);

% computing \Sigma matrix
var=2;
E=ABST.E;
xopt=ABST.xopt;
log_var=find(E.T==var);
Sigma=0;
for flcnt1=1:E.nsimple,
    nh=0;           % C-position counter for the right (C) terms 
    ng=0;           % C-position counter for the left (B) terms
    for flcnt2=1:size(E.X{E.simples(flcnt1)},2),
        nvar=E.X{E.simples(flcnt1)}(1,flcnt2);
        if nvar<0,
           X=decs2mats(xopt,E.L{log_var(abs(nvar))})';
        else
           X=decs2mats(xopt,E.L{log_var(nvar)});
        end
        L=E.B{E.simples(flcnt1)}(ng+1:ng+E.X{E.simples(flcnt1)}(2,flcnt2),:)';
        R=E.C{E.simples(flcnt1)}(nh+1:nh+E.X{E.simples(flcnt1)}(3,flcnt2),:);
        R=E.gqfm(flcnt1)*R;
        Sigma=Sigma+(L*X*R+R'*X'*L');
        nh=nh+E.X{E.simples(flcnt1)}(3,flcnt2);
        ng=ng+E.X{E.simples(flcnt1)}(2,flcnt2);
     end
end

% introduce the gain estimation terms
cz=ABST.c_z;
cf=ABST.c_f;
xg=ABST.xg;
Sigma=Sigma+cz'*cz-cf'*cf*xg;

% compute the main LMI
MainLMI=[P*A+A'*P,P*B;B'*P,zeros(m-n,m-n)]+Sigma;
