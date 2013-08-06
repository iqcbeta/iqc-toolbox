function [cost,xopt]=lmi_mincx_yalmip(obj)
% function [cost,xopt]=lmi_mincx_tbx(obj)
%
%  This program solves the LMI problem:
%
%       Minimize    c'x   subject to    L(x)  <  0
%       where L(x) is a matrix function. which is understood
%       as v'(L(x)+L(x)')v < 0 for all v in R^n.
%
%  The objective c'x (obj) has to be defined using the "abst"
%  environment "lmi". If the problem is feasible, a global cell
%  array ABSTSOLUTION, of same size as ABST.log, is created to
%  contain the optimal values of the decision variables.
%
% The default SDP solver is the solver coming with the yalmip
% IQC beta use YALMIP as the gateway to various SDP
% solvers. See YALMIP manual (http://users.isy.liu.se/johanl/yalmip/)
% for more details
%
% Last modified by cmj on 2013/4/29

global ABST
str={};

lmiparameter=ABST.lmiparameter;

% first, process the ABST log
if nargin == 0,
    E = lmi_extract;
    prob_feas = 1;
elseif nargin == 1,
    E = lmi_extract(obj);
    prob_feas = 0;
end

% symbolic names for interior types:
cst=1;
var=2;
lin=3;
lmi=4;

sft_v='(1e-9)';
% sht_v='0';
% sht_v='eps';

num_var_pos={};
F=[];

lmiparameter=ABST.lmiparameter;

%vrb=0 to suppress the output
vrb=1;
if ~isempty(ABST.lmiparameter)
    vrb = isempty(strfind(ABST.lmiparameter,'''verbose'',0'));
end

if ~vrb,
    disp_str(50,'lmi_mincx_yalmip')
else
    disp_str(51,'lmi_mincx_yalmip')
end

export_date=date;
str{1}=['%% Created by lmi_mincx_yalmip on ' export_date];
str{2}='\n%% Load coefficient matrices ...';
str{3}='load lmi_mincx_yalmip_exe';

%% initialize the YALMIP
str{4}='\n%% Initialize the YALMIP LMI';
str{5}='F = [];';
eval(str{5});

%% variable definite
if vrb,
    disp_str(52)
end
str{6}='\n%% Define multiplier variables ...';
sc=7; % string counter

vc=0;
kvar=find(E.T==var);
for k=1:length(kvar),
    ks   = num2str(k);
    varn = num2str(E.X{kvar(k)}(2,:));
    varm = num2str(E.X{kvar(k)}(3,:));
    switch ABST.log(kvar(k),5)
        case 1  % symmetric
            str{sc}=['X' ks '=sdpvar(' varn ',' varm ...
                ',''symmetric'');']; %#ok<*AGROW>
            eval(str{sc});
            sc=sc+1;
            num_var_pos=[num_var_pos;{['X',ks],'1',[]}];
        case 2  % rectangular
            str{sc}=['X' ks '=sdpvar(' varn ',' varm ',''full'');'];
            eval(str{sc});
            sc=sc+1;
            num_var_pos=[num_var_pos;{['X',ks],'2',[]}];
        case 3  % diagonal
            str{sc} = ['X' ks '=diag(sdpvar(' varn ',1));'];
            eval(str{sc});
            sc=sc+1;
            num_var_pos=[num_var_pos;{['X',ks],'3',[]}];
        case 4  % skew symmetric
            str{sc}=['X' ks '=sdpvar(' varn ',' varm ',''skew'');'];
            eval(str{sc});
            sc=sc+1;
            num_var_pos=[num_var_pos;{['X',ks],'4',[]}];
        case 5  % variable
            str{sc} = ['X' ks '=zeros(' varn ',' varm ');'];
            eval(str{sc});
            sc=sc+1;
            eval(['struct' ks '= E.L{' num2str(kvar(k)) '};']);
            M   = E.L{kvar(k)};
            Ma  = sort(M(:));
            Ms  = Ma(Ma>0);
            Msc = Ms(1);
            cnt = 1;
            for flcnt = 2: length(Ms)
                if Ms(flcnt) > Msc(cnt);
                    Msc(cnt+1) = Ms(flcnt);
                    cnt = cnt + 1;
                end
            end
            allvc=[];
            for flcnt = 1: length(Msc)
                Mscvar = num2str(Msc(flcnt));
                vc=vc+1;
                vcs = num2str(vc);
                str{sc}  = ['y' vcs '= sdpvar;'];
                eval(str{sc});
                sc=sc+1;
                str{sc} = ['X' ks '= X' ks '+((abs(struct' ks...
                    ')==' Mscvar ').*sign(struct' ks '))*y' vcs...
                    ';'];
                eval(str{sc});
                sc=sc+1;
                allvc=[allvc;['y',vcs]];
            end
            num_var_pos=[num_var_pos;{['X',ks],'5',...
                num2str(allvc)}];
        otherwise
            switch ABST.log(kvar(k),4),
                case 28, % subsref
                    rc=find(kvar==ABST.log(kvar(k),6));
                    p_rc=find(E.L{kvar(rc)}==E.L{kvar(k)});
                    str{sc} = ['X' ks '= X' num2str(rc)...
                        '(' num2str(p_rc(1)) ');' ];
                    eval(str{sc});
                    sc=sc+1;
                case 29, % subsasgn
                    rc=find(kvar==ABST.log(kvar(k),7));
                    p_rc=find(E.L{kvar(rc)}==E.L{kvar(k)});
                    str{sc} = ['X' ks '= X' num2str(rc)...
                        '(' num2str(p_rc(1)) ');' ];
                    eval(str{sc});
                    sc=sc+1;
                otherwise
                    disp_str(55)
            end
    end
end

%% now define all LMI's
lmi_log=find(E.T==lmi);
str{sc}='\n%% Define non-KYP lmis ...';
sc=sc+1;

for k=1:length(lmi_log),
    ks=num2str(k);
    kk=lmi_log(k);
    str{sc}=['LMI',ks,'=zeros(1);'];
    eval(str{sc});
    sc=sc+1;
    
    nh=0;
    ng=sum(E.X{kk}(3,:));
    eval(['CST',ks,'=[];']);
    for i=1:size(E.X{kk},2),
        is=num2str(i);
        vr=E.X{kk}(1,i);
        if isreal(vr)
            GG=E.C{kk}(ng+1:ng+E.X{kk}(2,i),:)';
            HH=E.C{kk}(nh+1:nh+E.X{kk}(3,i),:);
            eval(['l_',ks,'_',is,'=GG;'])
            eval(['r_',ks,'_',is,'=HH;'])
            if vr<0,
                vrs=num2str(-vr);
                str{sc}=['LMI',ks,'=LMI',ks,'+l_',ks,'_',is,'*X',vrs,...
                    '''*r_',ks,'_',is,'+(l_',ks,'_',is,'*X',vrs,...
                    '''*r_',ks,'_',is,')'';'];
                eval(str{sc});
                sc=sc+1;
            else
                vrs=num2str(vr);
                str{sc}=['LMI',ks,'=LMI',ks,'+l_',ks,'_',is,'*X',vrs,...
                    '*r_',ks,'_',is,'+(l_',ks,'_',is,'*X',vrs,...
                    '*r_',ks,'_',is,')'';'];
                eval(str{sc});
                sc=sc+1;
            end
        else
            GG=E.C{kk}(ng+1:ng+E.X{kk}(2,i),:)';
            HH=E.C{kk}(nh+1:nh+E.X{kk}(3,i),:);
            if isempty(eval(['CST',ks])),
                eval(['CST',ks,'=(GG*HH)+(GG*HH)'';']);
            else
                eval(['CST',ks,'=CST',ks,'+(GG*HH)+(GG*HH)'';']);
            end
        end
        nh=nh+E.X{kk}(3,i);
        ng=ng+E.X{kk}(2,i);
    end
    if ~isempty(eval(['CST',ks])),
        str{sc}=['LMI',ks,'=LMI',ks,'+CST',ks,';'];
        eval(str{sc})
        sc=sc+1;
    end
    eval(['n_lmi=size(LMI',ks,',1);'])
    str{sc} = ['F = F + [LMI',ks,'<=-',sft_v,...
        '*eye(',num2str(n_lmi),')];'];
    eval(str{sc});
    sc = sc+1;
end

%% now we define the "main" lmi
var_log=find(E.T==var);
if prob_feas==0,
    log_obj=double(obj);
    nh=0;
    ng=sum(E.X{log_obj}(3,:));
    count=0;
    COEFF=zeros(1,E.nlmivar+1);
    for i=1:size(E.X{log_obj},2),
        LM=E.C{log_obj}(ng+1:ng+E.X{log_obj}(2,i),:)';
        RM=E.C{log_obj}(nh+1:nh+E.X{log_obj}(3,i),:);
        
        if E.X{log_obj}(1,i)>0
            structureX=E.L{var_log(E.X{log_obj}(1,i))};
        elseif E.X{log_obj}(1,i)<0
            structureX=E.L{var_log(-E.X{log_obj}(1,i))}';
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
                for j=1:size(structureX,1)              % index of row of structureX
                    for s=1:size(structureX,2)          % index of column of structureX
                        nth_var=structureX(j,s);   % index of independent variables
                        if nth_var~=0,
                            if nth_var>E.nlmivar,
                                disp_str(66)
                            else
                                if nth_var>0
                                    COEFF(nth_var)=COEFF(nth_var)+conj(LM(j))*RM(s);
                                elseif nth_var<0
                                    COEFF(-nth_var)=COEFF(-nth_var)-conj(LM(j))*RM(s);
                                end
                            end
                        end
                    end
                end
            end
        end
        nh=nh+E.X{log_obj}(3,i);
        ng=ng+E.X{log_obj}(2,i);
    end
end

%% get variable
xnumber=[];
for i1=1:size(num_var_pos,1)
    eval(['[di_y,di_x]=size(',num_var_pos{i1,1},');']);
    switch num_var_pos{i1,2}
        case '1'
            for i2=1:di_x
                eval(['newxnumber=',num_var_pos{i1,1},'(',...
                    '1:',num2str(i2),',',num2str(i2),');']);
                newxnumber=reshape(newxnumber,1,size(newxnumber,1));
                xnumber=[xnumber,newxnumber];
            end
        case '2'
            for i2=1:di_y
                eval(['newxnumber=',num_var_pos{i1,1},'(',num2str(i2),...
                    ',:);']);
                xnumber=[xnumber,newxnumber];
            end
        case '3';
            for i2=1:di_x
                eval(['xnumber=[xnumber,',num_var_pos{i1,1},'(',...
                    num2str(i2),',',num2str(i2),')];']);
            end
        case '4';
            for i2=1:(di_x-1)
                eval(['xnumber=[xnumber,',num_var_pos{i1,1},'(',...
                    num2str(i2),',',num2str(i2+1),':',...
                    num2str(di_x),')];']);
            end
        case '5';
            for i2=1:size(num_var_pos{i1,3},1)
                eval(['newxnumber=num_var_pos{i1,3}(',...
                    num2str(i2),',:);']);
                eval(['xnumber=[xnumber,',newxnumber,'];']);
            end
    end
end
xnumber=reshape(xnumber,size(xnumber,2),1);

%% solving the system of LMI's
str{sc}='\n%% Solving the system of LMIs ...';
sc=sc+1;
if ~isempty(ABST.lmiparameter),
    str{sc}=['ops=sdpsettings(',ABST.lmiparameter,');'];
else
    str{sc}='ops=[];';
end
eval(str{sc});
sc=sc+1;

if prob_feas,
    if vrb,
        disp_str(56,num2str(length(xnumber)))
    end
    str{sc}='sol = solvesdp(F,[],ops);';
    eval(str{sc})
    sc=sc+1;
    
    checkset(F)
    [primalfeas,dualfeas]=checkset(F);
    switch sol.problem
        case 0 % success
            if vrb
                disp_str(67,'feasible')
            end
            xopt=double(xnumber);
            cost=-min(primalfeas);
        case 1
            if vrb
                disp_str(67,'infeasible')
            end
            xopt=[];
            cost=0;
        case {4,5}
            alleig=[];
            for i1=1:size(F)
                alleig=[alleig;eig(double(F{i1}+...
                    eval(sft_v)*eye(size(F{i1},1))))];
            end
            if min(alleig)>0
                if vrb
                    disp_str(67,'feasible');
                end
                xopt=double(xnumber);
                cost=-min(alleig);
            elseif min(alleig)<0 && -min(alleig) < str2double(sft_v)
                if vrb
                    disp_str(58,num2str(min(alleig)))
                end
                xopt=double(xnumber);
                cost=-min(alleig);
            else
                if vrb
                    disp_str(67,'infeasible')
                end
                xopt=[];
                cost=0;
            end
        otherwise
            if vrb
                disp_str(59,yalmiperror(sol.problem))
            end
            xopt=[];
            cost=0;
    end
else
    cmatrix=[COEFF(1:length(COEFF)-1),...
        zeros(1,length(xnumber)-length(COEFF)+1)];
    str{sc}='cmatrix_xnumber=cmatrix*xnumber;';
    eval(str{sc});
    sc=sc+1;
    if vrb,
        disp_str(56,num2str(length(xnumber)))
    end
    str{sc}='sol = solvesdp(F,cmatrix_xnumber,ops);';
    eval(str{sc})
    sc=sc+1;
    checkset(F)
    switch sol.problem
        case 0 % success
            if vrb
                disp_str(57,'Variable','feasible')
            end
            xopt=double(xnumber);
            str{sc}='cost=double(cmatrix_xnumber);';
            eval(str{sc})
            sc=sc+1;
        case 1
            if vrb
                disp_str(57,'Variable','infeasible')
            end
            xopt=[];
            cost=0;
        case {4,5}
            alleig=[];
            for i1=1:size(F)
                alleig=[alleig;eig(double(F{i1}+...
                    eval(sft_v)*eye(size(F{i1},1))))];
            end
            if min(alleig)>0
                if vrb
                    disp_str(57,'Variable','feasible');
                end
                xopt=double(xnumber);
                str{sc}='cost=double(cmatrix_xnumber);';
                eval(str{sc})
                sc=sc+1;
            elseif min(alleig)<0 && -min(alleig) < str2double(sft_v)
                if vrb
                    disp_str(58,num2str(min(alleig)))
                end
                xopt=double(xnumber);
                str{sc}='cost=double(cmatrix_xnumber);';
                eval(str{sc})
                sc=sc+1;
            else
                if vrb
                    disp_str(57,'Variable','infeasible')
                end
                xopt=[];
                cost=0;
            end
        otherwise
            if vrb
                disp_str(59,yalmiperror(sol.problem))
            end
            xopt=[];
            cost=0;
    end
    str{sc}='cost=cost+COEFF(length(COEFF));';
    eval(str{sc})
    sc=sc+1;
end

if ~vrb,
    save lmi_mincx_yalmip_exe lmiparameter CST*
    
    if ~prob_feas
        save lmi_mincx_yalmip_exe COEFF cmatrix -append
    end
    
    if exist('struct1','var')
        save lmi_mincx_yalmip_exe struct* -append
    end
    
    if exist('l_1_1','var')
        save lmi_mincx_yalmip_exe l_* -append
    end
    
    if exist('r_1_1','var')
        save lmi_mincx_yalmip_exe r_* -append
    end
    
    fid=fopen('lmi_mincx_yalmip_exe.m','wt');
    for i=1:sc-1,
        fprintf(fid,[str{i} '\n']);
    end
    fclose(fid);
end

ABST.xopt=xopt;
ABST.E=E;

if ~isempty(xopt)
    lmi_value;
end
