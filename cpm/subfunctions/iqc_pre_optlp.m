%%%%%%%%%%%%%%%%%%%%%%%
%%                   %%
%%   iqc_pre_optlp   %%
%%                   %%
%%%%%%%%%%%%%%%%%%%%%%%

function L=iqc_pre_optlp(f,z,options)
% function L=iqc_pre_optlp(f,z,options)
%
% Pre-processor for LP type optimizer of iqc 
% optimization problem
%
% f,y can be of type "inp" or "sgn" 
%
% Written by cykao@mit.edu on October 25, 1998

global ABST

% first, process the iqc abst log
E=iqc_extract(f,z);
E=iqc_reduce(E);
vrb=(ABST.lmiparameter(5)==0)|(ABST.lmiparameter(5)==777);

% symbolic names for interior types:
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

if nargin<2, error('two inputs required'); end
if (~isa(f,'abst'))|(~isa(z,'abst')),
   error('all inputs must be from the "abst" class')
end
if (E.T(double(f))~=sgn)&(E.T(double(f))~=inp),
   error('First argument must be input or signal')
end
if (E.T(double(z))~=sgn)&(E.T(double(z))~=inp),
   error('Second argument must be input or signal')
end

% ... Now, prepare all LMI's ... 
% All LMIs are represented to Sum(G*X*H)<0
% If the LMI is frequence depedent, L.lmi.A and L.lmi.B are not empty

if vrb,
   disp('  preparing the non-KYP LMIs ...')
end
lmi_log=find(E.T==lmi);
var_log=find(E.T==var);
for flcnt1=1:length(lmi_log),
    nth_lmi=lmi_log(flcnt1);
    [nn,mm]=size(E.AB{nth_lmi});         % dimensions
    if nn>0,
       L.lmi{flcnt1}.A=E.AB{nth_lmi}(1:nn,1:nn);
       L.lmi{flcnt1}.B=E.AB{nth_lmi}(1:nn,(nn+1):mm);
    else
       L.lmi{flcnt1}.A=[];
       L.lmi{flcnt1}.B=[];
    end
    nh=0;
    ng=sum(E.X{nth_lmi}(3,:));
    count=0;
    for flcnt2=1:size(E.X{nth_lmi},2),
        if E.X{nth_lmi}(1,flcnt2)>0
           structureX=E.L{var_log(E.X{nth_lmi}(1,flcnt2))};
        elseif E.X{nth_lmi}(1,flcnt2)<0
           structureX=E.L{var_log(-E.X{nth_lmi}(1,flcnt2))}';
        end
        LM=E.C{nth_lmi}(ng+1:ng+E.X{nth_lmi}(2,flcnt2),:)';
        RM=E.C{nth_lmi}(nh+1:nh+E.X{nth_lmi}(3,flcnt2),:);

        L.lmi{flcnt1}.G{count+1} = LM;
        L.lmi{flcnt1}.X{count+1} = structureX;
        L.lmi{flcnt1}.H{count+1} = RM;
        count=count+1;

        nh=nh+E.X{nth_lmi}(3,flcnt2);
        ng=ng+E.X{nth_lmi}(2,flcnt2);
    end
    L.lmi{flcnt1}.size=size(L.lmi{flcnt1}.G{1},1);
    L.lmi{flcnt1}.mlmi=0;
end

if vrb,
   disp('  preparing the KYP LMIs ...')
end

mainlmi=length(lmi_log)+1;

if E.nstates>0,
   L.lmi{mainlmi}.A=E.ab(1:E.nstates,1:E.nstates);
   L.lmi{mainlmi}.B=E.ab(1:E.nstates,(E.nstates+1):size(E.ab,2));
else
   L.lmi{mainlmi}.A=[];
   L.lmi{mainlmi}.B=[];
end
L.lmi{mainlmi}.size=size(E.ab,2);
L.lmi{mainlmi}.mlmi=1;
count=0;
for flcnt1=1:E.nsimple,
    nh=0;           % C-position counter for the right (C) terms
    ng=0;           % C-position counter for the left (B) terms
    for flcnt2=1:size(E.X{E.simples(flcnt1)},2),
        nth_var=E.X{E.simples(flcnt1)}(1,flcnt2);
        if nth_var>0,
           structureX=E.L{var_log(nth_var)};
        elseif nth_var<0,
           structureX=E.L{var_log(-nth_var)}';
        end
        LM=E.B{E.simples(flcnt1)}(ng+1:ng+E.X{E.simples(flcnt1)}(2,flcnt2),:)';
        RM=E.C{E.simples(flcnt1)}(nh+1:nh+E.X{E.simples(flcnt1)}(3,flcnt2),:);
        RM=E.gqfm(flcnt1)*RM;

        L.lmi{mainlmi}.G{count+1}=LM;
        L.lmi{mainlmi}.X{count+1}=structureX;
        L.lmi{mainlmi}.H{count+1}=RM;
        count=count+1;

        nh=nh+E.X{E.simples(flcnt1)}(3,flcnt2);
        ng=ng+E.X{E.simples(flcnt1)}(2,flcnt2);
    end
end

% introduce the gain estimation terms
switch options,
  case 1,    % L2 gain estimation
    c_z=E.C{double(z)};
    c_f=E.C{double(f)};
    L.lmi{mainlmi}.G{count+1}=c_f';
    L.lmi{mainlmi}.X{count+1}=sum(E.nlmivar)+1;
    L.lmi{mainlmi}.H{count+1}=-c_f;
    L.lmi{mainlmi}.G{count+2}=eye(length(c_z));
    L.lmi{mainlmi}.X{count+2}=[];
    L.lmi{mainlmi}.H{count+2}=eye(length(c_z));
    L.lmi{mainlmi}.C=c_z'*c_z;
    L.ndecvar=sum(E.nlmivar)+1;
    L.obj=[zeros(1,sum(E.nlmivar)),1];

  case 2,    % L2 --> L_infinity estimation
    c_f=E.C{double(f)};
    c_z=E.C{double(z)};
    ns=size(E.ab,1);
    cmatrix_z=c_z(:,1:ns);
    dmatrix_z=c_z(:,ns+1:size(c_z,2));
    if any(any(dmatrix_z))
       error('L2 --> L_infinity estimation : need output signal to be differentiable!')
    else
       c_z_dot=cmatrix_z*E.ab;
       L.lmi{mainlmi}.G{count+1}=c_f';
       L.lmi{mainlmi}.X{count+1}=sum(E.nlmivar)+1;
       L.lmi{mainlmi}.H{count+1}=-c_f;
       L.lmi{mainlmi}.G{count+2}=eye(length(c_z));
       L.lmi{mainlmi}.X{count+2}=[];
       L.lmi{mainlmi}.H{count+2}=eye(length(c_z_dot));
       L.lmi{mainlmi}.C = c_z'*c_z_dot + c_z_dot'*c_z;
       L.ndecvar=sum(E.nlmivar)+1;
       L.obj=[zeros(1,sum(E.nlmivar)),1];
    end
end

global ABST
ABST.E=E;
ABST.L=L;
