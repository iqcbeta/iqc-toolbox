function L=lmi_lp_pre1(E,obj)
% function L=lmi_lp_pre1(E,obj)
%
% Preparation for LP based LMI solver. Version 1
% 
% written by cykao@mit.edu, last modified: Feb. 09 1999
%                           last modified: Apr. 29 1999

% safety checks
global ABST
A=ABST;
if ~isfield(A,'name'), 
   error('"abst" environment not defined')
end
if ~strcmp(A.name,'lmi') & ~strcmp(A.name,'fdlmi'),
   error('This is not a "lmi" or a "fdlmi" environment')
end

if nargin==1, prob_feas=1; else prob_feas=0; end

% symbolic names for interior types:
cst=1;
var=2;
lin=3;
lmi=4;

% Constructing the "L" field
global L
lmi_log=find(E.T==lmi);
tmp1=[lmi_log,A.log(lmi_log,2)];    % sort the lmi_log by the size of LMIs
tmp1=sortrows(tmp1,2);
lmi_log([length(lmi_log):-1:1],1)=tmp1(:,1);

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
    CST=[];
    for flcnt2=1:size(E.X{nth_lmi},2),

        LM=E.C{nth_lmi}(ng+1:ng+E.X{nth_lmi}(2,flcnt2),:)';
        RM=E.C{nth_lmi}(nh+1:nh+E.X{nth_lmi}(3,flcnt2),:);

        if E.X{nth_lmi}(1,flcnt2)>0
           structureX=E.L{var_log(E.X{nth_lmi}(1,flcnt2))};
           L.lmi{flcnt1}.G{count+1} = LM;
           L.lmi{flcnt1}.X{count+1} = structureX;
           L.lmi{flcnt1}.H{count+1} = RM;
           count=count+1;
        elseif E.X{nth_lmi}(1,flcnt2)<0
           structureX=E.L{var_log(-E.X{nth_lmi}(1,flcnt2))}';
           L.lmi{flcnt1}.G{count+1} = LM;
           L.lmi{flcnt1}.X{count+1} = structureX;
           L.lmi{flcnt1}.H{count+1} = RM;
           count=count+1;
        else
           if isempty(CST)
              CST=(LM*RM);
           else
              CST=CST+(LM*RM);
           end
        end
        nh=nh+E.X{nth_lmi}(3,flcnt2);
        ng=ng+E.X{nth_lmi}(2,flcnt2);
    end
    L.lmi{flcnt1}.size=size(L.lmi{flcnt1}.G{1},1);
    if ~isempty(CST)
        L.lmi{flcnt1}.C=CST;
    else
        L.lmi{flcnt1}.C=zeros(L.lmi{flcnt1}.size);
    end
    L.lmi{flcnt1}.G{count+1} = eye(size(L.lmi{flcnt1}.G{1},1));
    L.lmi{flcnt1}.X{count+1} = [];
    L.lmi{flcnt1}.H{count+1} = eye(size(L.lmi{flcnt1}.G{1},1));
end

L.ndecvar = E.nlmivar;

% prepare the cost function
if prob_feas==0,
   log_obj=double(obj);
   nh=0;
   ng=sum(E.X{log_obj}(3,:));
   count=0;
   COEFF=zeros(1,E.nlmivar);         
   
   for flcnt2=1:size(E.X{log_obj},2),
       LM=E.C{log_obj}(ng+1:ng+E.X{log_obj}(2,flcnt2),:)';
       RM=E.C{log_obj}(nh+1:nh+E.X{log_obj}(3,flcnt2),:);
       
       if E.X{log_obj}(1,flcnt2)>0
          structureX=E.L{var_log(E.X{log_obj}(1,flcnt2))};
       elseif E.X{log_obj}(1,flcnt2)<0
          structureX=E.L{var_log(-E.X{log_obj}(1,flcnt2))}';
       else
          structureX=[];
          CSTD=E.B{E.X{log_obj}(1,flcnt2)*(-sqrt(-1))};
          COEFF(1,E.nlmivar+1)=COEFF(1,E.nlmivar+1)+LM*CSTD*RM;
       end

       if ~isempty(structureX)
          mX=size(structureX,1);
          nX=size(structureX,2);
          if (mX==1 & nX==1)
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
                           error('Bad index of independent variables!!')
                        else
                           if nth_var>0
                              COEFF(nth_var)=COEFF(nth_var)+LM(flcnt3)*RM(flcnt4);
                           elseif nth_var<0
                              COEFF(-nth_var)=COEFF(-nth_var)-LM(flcnt3)*RM(flcnt4);
                           end
                        end
                     end
                 end
             end
          end
       end

       nh=nh+E.X{log_obj}(3,flcnt2);
       ng=ng+E.X{log_obj}(2,flcnt2);
   end
   L.obj=COEFF;
else
   L.obj=[];
end
