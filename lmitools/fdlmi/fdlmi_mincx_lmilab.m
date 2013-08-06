function [cost,xopt]=fdlmi_mincx_lmilab(obj)
% function [cost,xopt]=fdlmi_mincx_YALMIP(obj)
%
%  This program solves the LMI problem:
%
%       Minimize    c'x   subject to    L(w,x)  <  0
%       where L(w,x) is a matrix function, which can be frequency dependent.
%       If L is frequency dependent, then L(w,x) < 0 should be understood
%       as v'Re(L(w,x))v < 0 for all w, and for all v in R^n (w is the frequency)
%
%  The objective c'x (obj) has to be defined using the "abst"
%  environment "fdlmi". If the problem is feasible, a global cell
%  array ABSTSOLUTION, of same size as ABST.log, is created to
%  contain the optimal values of the decision variables.
%
%  This function uses the LMI Control Toolbox by The MathWorks, Inc.
%
%  Written by cykao@mit.edu    Mar. 15 1999
%             last modified    Apr. 29 1999
%             last modified    Jul. 13 1999
% Last modified by cmj on 2013/4/27

global ABST
str={};

ABST.systemtype = 'continuous';

lmiparameter=ABST.lmiparameter;

% first, process the ABST log
if nargin == 0,
    E = fdlmi_extract;
    prob_feas = 1;
elseif nargin == 1,
    E = fdlmi_extract(obj);
    prob_feas = 0;
end

% symbolic names for interior types:
cst=1;
var=2;
lin=3;
lmi=4;

% vrb=0 to suppress the output
vrb=(lmiparameter(5)==0)||(lmiparameter(5)==777);

if vrb,
    if lmiparameter(5)>100,
        disp_str(50,'fdlmi_mincx_lmilab')
    else
        disp_str(51,'fdlmi_mincx_lmilab')
    end
end

export_date=date;
str{1}=['%% Created by fdlmi_mincx_lmilab ' export_date];
str{2}='\n%% Load coefficient matrices ...';
str{3}='load fdlmi_mincx_lmilab_exe';

%% initialize the LMI Control Toolbox
str{4}='\n%% Initialize the LMI Lab';
str{5}='setlmis([]);';
eval(str{5})

%% now define all variables
if vrb,
    disp_str(52)
end

kvar=find(E.T==var);
str{6}='\n%% Define multiplier variables ...';
sc=7;              % string counter
for k=1:length(kvar),
    ks=num2str(k);
    eval(['struct',ks,'=E.L{' num2str(kvar(k)),'};'])
    str{sc}=['x',ks,'=lmivar(3,struct' ks ');']; %#ok<*AGROW>
    eval(str{sc});
    sc=sc+1;
end

%% now we define the "non-KYP" lmi
lmi_log=find(E.T==lmi);
str{sc}='\n%% Define non-KYP lmis ...';
sc=sc+1;
freq_dep_lmi_flag=0;

for k=1:length(lmi_log),
    ks=num2str(k);
    kk=lmi_log(k);
    [nn,mm]=size(E.AB{kk});    % dimensions
    nns=num2str(nn);
    vs=mm-nn;                       % LMI size
    vss=num2str(vs);
    if nn>0,
        freq_dep_lmi_flag=1;
        eval(['state',ks,'=E.AB{kk};']);
        str{sc}=['p',ks,'=lmivar(1,[',nns,' 1]);'];
        eval(str{sc});
        sc=sc+1;
        if strcmp(ABST.systemtype,'continuous'),
            str{sc}=['lmiterm([',ks,' 1 1 p',ks,...
                '],[eye(',nns,');zeros(',vss,',',nns,...
                ')],state',ks,',''s'');'];
            eval(str{sc});
            sc=sc+1;
        elseif strcmp(ABST.systemtype,'discrete'),
            str{sc}=['lmiterm([',ks,' 1 1 p',ks,...
                '],state',ks,''',state',ks,');'];
            eval(str{sc});
            sc=sc+1;
            str{sc}=['lmiterm([',ks,' 1 1 p',ks,...
                '],[-eye(',nns,') zeros(',nns,',',vss,...
                ')]'',[eye(',nns,') zeros(',nns,',',vss,...
                ')]);'];
            eval(str{sc});
            sc=sc+1;
        else
            disp_str(55)
        end
    end
    nh=0;
    ng=sum(E.X{kk}(3,:));
    
    eval(['CST',ks,'=[];']);
    for i1=1:size(E.X{kk},2),
        nth_var=E.X{kk}(1,i1);
        if isreal(nth_var)
            is=num2str(i1);
            vrs=num2str(nth_var);   % variable number
            GG=E.C{kk}(ng+1:ng+E.X{kk}(2,i1),:)';
            HH=E.C{kk}(nh+1:nh+E.X{kk}(3,i1),:);
            eval(['l_',ks,'_',is,'=GG;'])
            eval(['r_',ks,'_',is,'=HH;'])
            str{sc}=['lmiterm([',ks,' 1 1 ',vrs,'],l_',ks,...
                '_',is,',r_',ks,'_',is,',''s'');'];
            eval(str{sc});
            sc=sc+1;
        else
            GG=E.C{kk}(ng+1:ng+E.X{kk}(2,i1),:)';
            HH=E.C{kk}(nh+1:nh+E.X{kk}(3,i1),:);
            if isempty(eval(['CST',ks]))
                eval(['CST',ks,'=(GG*HH)+(GG*HH)'';']);
            else
                eval(['CST',ks,'=CST',ks,'+(GG*HH)+(GG*HH)'';']);
            end
        end
        nh=nh+E.X{kk}(3,i1);
        ng=ng+E.X{kk}(2,i1);
    end
    if ~isempty(eval(['CST',ks])),
        str{sc}=['lmiterm([',ks,' 1 1 0],CST',ks,');'];
        eval(str{sc})
        sc=sc+1;
    end
end

%% now we define the "main" lmi
var_log=find(E.T==var);
if prob_feas==0,
    log_obj=double(obj);
    nh=0;
    ng=sum(E.X{log_obj}(3,:));
    count=0;
    COEFF=zeros(1,E.nlmivar+1);
    for i1=1:size(E.X{log_obj},2),
        LM=E.C{log_obj}(ng+1:ng+E.X{log_obj}(2,i1),:)';
        RM=E.C{log_obj}(nh+1:nh+E.X{log_obj}(3,i1),:);
        
        if E.X{log_obj}(1,i1)>0
            structureX=E.L{var_log(E.X{log_obj}(1,i1))};
        elseif E.X{log_obj}(1,i1)<0
            structureX=E.L{var_log(-E.X{log_obj}(1,i1))}';
        else
            structureX=[];
            COEFF(1,E.nlmivar+1)=COEFF(1,E.nlmivar+1)+LM*RM;
        end
        
        if ~isempty(structureX)
            mX=size(structureX,1);
            nX=size(structureX,2);
            if (mX==1 && nX==1)
                nth_var=structureX(1,1);
                if nth_var>0
                    COEFF(nth_var)=COEFF(nth_var)+LM*RM;
                elseif nth_var<0
                    COEFF(-nth_var)=COEFF(-nth_var)-LM*RM;
                end
            else
                for flcnt3=1:size(structureX,1)              % index of row of structureX
                    for flcnt4=1:size(structureX,2)          % index of column of structureX
                        nth_var=structureX(flcnt3,flcnt4);   % index of independent variables
                        if nth_var~=0,
                            if nth_var>E.nlmivar,
                                disp_str(66)
                            else
                                if nth_var>0
                                    COEFF(nth_var)=COEFF(nth_var)+...
                                        conj(LM(flcnt3))*RM(flcnt4);
                                elseif nth_var<0
                                    COEFF(-nth_var)=COEFF(-nth_var)-conj(LM(flcnt3))*RM(flcnt4);
                                end
                            end
                        end
                    end
                end
            end
        end
        nh=nh+E.X{log_obj}(3,i1);
        ng=ng+E.X{log_obj}(2,i1);
    end
end

%% solving the system of LMI's
str{sc}='\n%% Solving the system of LMIs ...';
sc=sc+1;
str{sc}='lmi=getlmis;';
eval(str{sc})
sc=sc+1;
str{sc}='ndecvar=decnbr(lmi);';
eval(str{sc})
sc=sc+1;

cost=0;
xopt=[];

if prob_feas
    if vrb
        disp_str(56,num2str(ndecvar))
    end
    cost=0;
    str{sc}='[cost,xopt] = feasp(lmi,lmiparameter);';
    eval(str{sc})
    sc=sc+1;
    if cost<0
        if vrb,
            disp_str(67,'feasible');
        end
    else
        if vrb,
            disp_str(67,'infeasible');
        end
    end
else
    cmatrix=[COEFF(1:length(COEFF)-1),zeros(1,ndecvar-length(COEFF)+1)];
    if vrb,
        disp_str(56,num2str(ndecvar))
    end
    str{sc}='[cost,xopt]=mincx(lmi,cmatrix,lmiparameter);';
    eval(str{sc})
    sc=sc+1;
    if ~isempty(cost)
        if vrb,
            disp_str(57,'Variable','feasible')
        end
        str{sc}='cost=cost+COEFF(length(COEFF));';
        eval(str{sc})
        sc=sc+1;
    else
        if vrb,
            disp_str(57,'Variable','infeasible')
        end
        cost=0;
    end
    
end

%%
if lmiparameter(5)>100
    save fdlmi_mincx_lmilab_exe lmiparameter CST*
    
    if ~prob_feas
        save fdlmi_mincx_lmilab_exe COEFF cmatrix -append
    end
    
    if freq_dep_lmi_flag
        save fdlmi_mincx_lmilab_exe state* -append
    end
    
    if exist('struct1','var')
        save fdlmi_mincx_lmilab_exe struct* -append
    end
    
    if exist('l_1_1','var')
        save fdlmi_mincx_lmilab_exe l_* -append
    end
    
    if exist('r_1_1','var')
        save fdlmi_mincx_lmilab_exe r_* -append
    end
    
    fid=fopen('fdlmi_mincx_lmilab_exe.m','wt');
    for i1=1:sc-1,
        fprintf(fid,[str{i1} '\n']);
    end
    fclose(fid);
end

ABST.xopt=xopt;
ABST.E=E;

if ~isempty(xopt),
    fdlmi_value;
end
