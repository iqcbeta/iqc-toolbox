function [cost,xopt]=fdlmi_mincx_tbx_d(obj)
% function [cost,xopt]=fdlmi_mincx_tbx(obj)
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
%  Written by cykao@mit.edu  Mar. 15 1999
%             last modified  Apr. 29 1999
%             last modified  Jul. 13 1999


% first, process the ABST log
% ---------------------------
if nargin==0,
   E=fdlmi_extract;
   prob_feas=1;
elseif nargin==1,
   E=fdlmi_extract(obj);
   prob_feas=0;
end

% symbolic names for interior types:
cst=1;
var=2;
lin=3;
lmi=4;

% vrb=0 to suppress the output
vrb=(E.options(5)==0)|(E.options(5)==777);

if vrb,
   if E.options(5)>100,
      disp('iqc_gain_tbx  (exporting the LMI script to lmi_exe.m) ...')
   else
      disp('iqc_gain_tbx ...')
   end
end

export_date=date;
str{1}=['%% Created by iqc_gain_tbx on ' export_date];
str{2}='\n%% Load coefficient matrices ...';
str{3}=['load lmi_exe'];

% initialize the LMI Control Toolbox
str{4}='\n%% Initialize the LMI Lab';
str{5}='setlmis([]);';   
eval(str{5})

% now define all variables
if vrb,
   disp('  defining the original variables ...')
end

kvar=find(E.T==var);
str{6}='\n%% Define multiplier variables ...';
sc=7;              % string counter
for flcnt1=1:length(kvar),
    ks=num2str(flcnt1);
    eval(['struct' ks '=E.L{' num2str(kvar(flcnt1)) '};'])
    str{sc}=['x' ks '=lmivar(3,struct' ks ');'];
    eval(str{sc});
    sc=sc+1;
end

lmi_log=find(E.T==lmi);
str{sc}='\n%% Define non-KYP lmis ...';
sc=sc+1;
freq_dep_lmi_flag=0;
for flcnt1=1:length(lmi_log),
    nth_lmi=lmi_log(flcnt1);
    [mAB,nAB]=size(E.AB{nth_lmi});    % dimensions
    ks=num2str(flcnt1);
    nns=num2str(mAB);
    vs=nAB-mAB;                       % LMI size
    vss=num2str(vs);
    if mAB>0,
       freq_dep_lmi_flag=1;
       eval(['state' ks '=E.AB{nth_lmi};']);
       str{sc}=['P' ks '=lmivar(1,[' nns ' 1]);'];
       eval(str{sc});
       sc=sc+1;
       str{sc}=['lmiterm([',ks,' 1 1 P',ks,'], state',ks,''', state',ks,');'];
       eval(str{sc});
       sc=sc+1;
       str{sc}=['lmiterm([',ks,' 1 1 P',ks,'], [-eye(',nns,');zeros(',vss,',',nns,')], [eye(',nns,'),zeros(',nns,',',vss,')] );'];
       eval(str{sc});
       sc=sc+1;
    end
    nh=0;
    ng=sum(E.X{nth_lmi}(3,:));
    CST=[];
    for flcnt2=1:size(E.X{nth_lmi},2),
        nth_var=E.X{nth_lmi}(1,flcnt2);
        if isreal(nth_var)
           is=num2str(flcnt2);
           nvs=num2str(nth_var);   % variable number
           GG=E.C{nth_lmi}(ng+1:ng+E.X{nth_lmi}(2,flcnt2),:)';
           HH=E.C{nth_lmi}(nh+1:nh+E.X{nth_lmi}(3,flcnt2),:);
           eval(['l_' ks '_' is '=GG;'])
           eval(['r_' ks '_' is '=HH;'])
           str{sc}=['lmiterm([' ks ' 1 1 ' nvs '],l_' ks ...
                     '_' is ',r_' ks '_' is ',''s'');'];
           eval(str{sc});
           sc=sc+1;
        else
           % nth_var=nth_var*(-sqrt(-1));
           GG=E.C{nth_lmi}(ng+1:ng+E.X{nth_lmi}(2,flcnt2),:)';
           HH=E.C{nth_lmi}(nh+1:nh+E.X{nth_lmi}(3,flcnt2),:);
           if isempty(CST)
              CST=(GG*HH)+(GG*HH)';
           else
              CST=CST+(GG*HH)+(GG*HH)';
           end
        end
        nh=nh+E.X{nth_lmi}(3,flcnt2);
        ng=ng+E.X{nth_lmi}(2,flcnt2);
    end
    if ~isempty(CST)
        str{sc}=['lmiterm([' ks ' 1 1 0],CST);'];
        eval(str{sc})
        sc=sc+1;
    end
end

var_log=find(E.T==var);
if prob_feas==0,
   log_obj=double(obj);
   nh=0;
   ng=sum(E.X{log_obj}(3,:));
   count=0;
   COEFF=zeros(1,E.nlmivar+1);         
   for flcnt2=1:size(E.X{log_obj},2),
       LM=E.C{log_obj}(ng+1:ng+E.X{log_obj}(2,flcnt2),:)';
       RM=E.C{log_obj}(nh+1:nh+E.X{log_obj}(3,flcnt2),:);
       
       if E.X{log_obj}(1,flcnt2)>0
          structureX=E.L{var_log(E.X{log_obj}(1,flcnt2))};
       elseif E.X{log_obj}(1,flcnt2)<0
          structureX=E.L{var_log(-E.X{log_obj}(1,flcnt2))}';
       else
          structureX=[];
          COEFF(1,E.nlmivar+1)=COEFF(1,E.nlmivar+1)+LM*RM;
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
                              COEFF(nth_var)=COEFF(nth_var)+conj(LM(flcnt3))*RM(flcnt4);
                           elseif nth_var<0
                              COEFF(-nth_var)=COEFF(-nth_var)-conj(LM(flcnt3))*RM(flcnt4);
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
end
 
% solving the system of LMI's 
str{sc}='\n%% Solving the system of LMIs ...';
sc=sc+1;
str{sc}='lmi=getlmis;';
eval(str{sc})
str{sc+1}='ndecvar=decnbr(lmi);';
eval(str{sc+1})
sc=sc+1;

options=E.options;
if prob_feas
   str{sc}='[cost,xopt] = feasp(lmi,options);';
   sc=sc+1;
   if vrb,
      disp(['  Solving with ' num2str(decnbr(lmi)) ' decision variables ...'])
   end
   eval(str{sc-1})
   disp(' ')
   if cost<0
      if vrb,
         disp('The problem is feasible.')
      end
   else
      if vrb,
         disp('The problem is infeasible.')
      end
   end
else
   cmatrix=[COEFF(1:length(COEFF)-1),zeros(1,ndecvar-length(COEFF)+1)];
   str{sc}='[cost,xopt]=mincx(lmi,cmatrix,options);';
   sc=sc+1;
   str{sc}='cost=cost+COEFF(length(COEFF));';
   if vrb,
      disp(['  Solving with ' num2str(decnbr(lmi)) ' decision variables ...'])
   end
   eval(str{sc-1})
   eval(str{sc})
end

if E.options(5)>100,
   save lmi_exe options cmatrix COEFF CST

   if freq_dep_lmi_flag
      save lmi_exe state* -append
   end

   if exist('struct1')
      save lmi_exe struct* -append
   end

   if exist('l_1_1')
      save lmi_exe l_* -append
   end

   if exist('r_1_1')
      save lmi_exe r_* -append
   end
   
   fid=fopen('lmi_exe.m','wt');
   for i=1:sc,
      fprintf(fid,[str{i} '\n']);
   end
   fclose(fid);
end


global ABST
ABST.E=E;
ABST.xopt=xopt;
if ~isempty(xopt),
   fdlmi_value;
end
