function [cost,xopt,E,L]=fdlmi_mincx_lp2(obj,method)
% function [cost,xopt]=fdlmi_mincx_lp2(obj,method)
%
%  This program solves the LMI problem:
%
%       Minimize    c'x   subject to    L(w,x)  <  0
%       where L(w,x) is a matrix function, which can be frequency dependent
%       If L are frequency dependent, then L(w,x) < 0 should be understood
%       as v'Re(L(w,x))v < 0 for all w, and for all v in R^n (w is the frequency)
%
%  by a cutting plane algorithm developed by Megretski's group at
%  LIDS M.I.T.
%
% Input arguments:
%   'obj' is the cost function (c'x) that will be minimized
%    note :'obj' must be class "abst" (var or lin)
%
%   'method' specifies which LP solver will be used in solving
%    the optimization problem:
%
%      'm' : using MATLAB LP solver lp.m
%      'a' : using the structured linear program solver written by
%            Prof. A Megretski (ameg@mit.edu)
%       default : method='a'
%
% Output arguments:
%    'cost' : optimal objective
%    'xopt' : optimal decision variables that acheieve the optimal objective
%
% written by cykao@mit.edu  Feb. 09 1999
%            last modified  Apr. 30 1999

% initialization
% --------------
if nargin==0,
    E=fdlmi_extract;
    lmi_lp_pre1(E);
    method='a';
elseif nargin==1,
    E=fdlmi_extract(obj);
    lmi_lp_pre1(E,obj);
    method='a';
else
    E=fdlmi_extract(obj);
    lmi_lp_pre1(E,obj);
end

global LP_A LP_B L %#ok<*REDEF>
global ndecv_e ndecv_o
global alpha
global error_tol1 error_tol2
global no_improve
global more_constraint
global re_solve_LP
global feasible

max_iteration=250;              % maximum iteration
LP_C=[zeros(1,E.nlmivar),1];    % set LP_C, LP_A, LP_B for feasibility problem
LP_A=[zeros(1,E.nlmivar),-1];
LP_B=1;
ndecv_e=length(LP_C);           % total number of decision variables
ndecv_o=ndecv_e-1;              % number of decision variables from LMI problem
y_ubound=[];                    % upper bound of max eigenvalue
alpha=0.382;
error_tol1=1e-5;
error_tol2=1e-5;                % error torlerance of (sub-optimal)-(true optimal).
no_improve=0;

% find out the non-frequence depedent LMIs
% and frequence depedent LMIs
% ----------------------------------------
fdeplmi=[];
nfdeplmi=[];
phc_A=[];
phc_B=[];
for flcnt1=1:length(L.lmi) %#ok<*NODEF>
    if isempty(L.lmi{flcnt1}.A)
        if L.lmi{flcnt1}.size==1,  % collect original polyhedral constraints
            [coeff,constant]=hyperplane(L,flcnt1,1,ndecv_o);
            phc_A=[phc_A;[coeff,-1]]; %#ok<*AGROW>
            phc_B=[phc_B;-constant];
            L.lmi{flcnt1}.done=1; %#ok<*STRNU>
        else
            nfdeplmi=[nfdeplmi;flcnt1];
        end
    else
        fdeplmi=[fdeplmi;flcnt1];
    end
end

LP_A=[LP_A;phc_A];
LP_B=[LP_B;phc_B];

% ------------------- %
%     Main Loop       %
% ------------------- %

done=0;          % flag : set to 1 when optimal solution found
counter=1;       % counter for how many loops the program has been running
feasible=0;      % feasibility flag; feasible=1 if feasible point found.
re_solve_LP=1;   % flag : set to 1 when there is a need to re-solve LP

disp(' ')
disp('Now, solving LMIs ......')
disp(' ')
%%             1         2         3         4
%% ## 1234567890123456789012345678901234567890
disp(' Iteration    :     y upper bound      ')
disp(' ')

while ~done && counter<=max_iteration,
    
    % ------------------------ %
    %  INITILAIZATION :        %
    %  randomly generating     %
    %  initial LP constraints  %
    % ------------------------ %
    
    if counter==1,
        y_ubound=initialize(nfdeplmi,fdeplmi);
    end
    
    % ---------------- %
    %  solve LP ...    %
    % ---------------- %
    
    if re_solve_LP==1,
        [lambda,y]=LP_solve(LP_C,method);
        if y > 0  %%% infeasibility estabilished !! %%%
            gain=[];
            xopt=[];
            disp('The problem is infeasible')
            return
        end
    end
    
    % compute \bar(y)
    y_bar=alpha*y+(1-alpha)*y_ubound;
    
    re_solve_LP=0;
    more_constraint=0;  % A flag that indicates whether there are more constraints formed
    % '0': no more constraints formed; '1': more constraints
    
    % ------------------------------------------------------------------- %
    %   check if the non-frequence depedent LMIs are satisfied; if not,   %
    %   computing the corresponding constraints on LP.                    %
    % ------------------------------------------------------------------- %
    for flcnt1=1:length(nfdeplmi)
        nth_lmi=nfdeplmi(flcnt1);
        MX=lambda_eval(lambda,nth_lmi);
        
        if feasible,
            MF=L.lmi{nth_lmi}.MF;
            MX_r=MX+y_bar*MF;                          % S(lambda)+y*S(lambda_0)<0 ?
        else
            L.lmi{nth_lmi}.MF=MX;                      % replace the S(lambda_0) by current value
            MX_r=MX-y_bar*eye(L.lmi{nth_lmi}.size);    % S(lambda)-y*I<0 ?
        end
        
        max_MX_eig=V_NFDLMI(MX_r,nth_lmi);
        
        if more_constraint==1,
            if feasible,
                yh_u=y_ubound;
                yh_l=y_bar;
                dist=abs(y_ubound-y_bar);
                beta=0.5;
                while abs(yh_u-yh_l)>=0.05*dist
                    yh_m=beta*yh_l+(1-beta)*yh_u;
                    MX_r=MX+yh_m*MF;
                    max_MX_eig=V_NFDLMI(MX_r,nth_lmi);
                    if more_constraint==0,
                        yh_u=yh_m;
                    elseif more_constraint==1,
                        yh_l=yh_m;
                    end
                end
                
                % re-compute Q,F,R by the new testing point
                y_bar=yh_u;
                more_constraint=0;
            else
                y_bar=y_bar+(max_MX_eig+1e-6);            % move y_bar back!
                more_constraint=0;
            end
        end
    end
    
    % ------------------------------------------------------------- %
    %    check if the frequency depedent LMIs are satisfied;        %
    %    if not, computing the corresponding constraints on LP.     %
    % ------------------------------------------------------------- %
    
    for flcnt1=1:length(fdeplmi)
        nth_lmi=fdeplmi(flcnt1);
        
        % ----------------------------------------------------- %
        %   Step 1: evaluating and forming Hamiltonian matrix   %
        % ----------------------------------------------------- %
        MX=lambda_eval(lambda,nth_lmi);
        
        A=L.lmi{nth_lmi}.A;                          % extract A, B, Q, F, R matrices
        B=L.lmi{nth_lmi}.B;
        nA=size(A,1);
        mA=size(A,2);
        mB=size(B,2);
        Q_o=MX(1:nA,1:mA);
        F_o=MX(1:nA,mA+1:mA+mB);
        R_o=MX(nA+1:nA+mB,mA+1:mA+mB);
        
        if feasible,
            Q_r=L.lmi{nth_lmi}.MF(1:nA,1:mA);
            F_r=L.lmi{nth_lmi}.MF(1:nA,mA+1:mA+mB);
            R_r=L.lmi{nth_lmi}.MF(nA+1:nA+mB,mA+1:mA+mB);
            Q=Q_o+y_bar*Q_r;                          % [Q,F;F',R]=S(lambda)+y*S(lambda_0)
            F=F_o+y_bar*F_r;
            R=R_o+y_bar*R_r;
        else
            L.lmi{nth_lmi}.MF=MX;                     % replace S(lambda_0) by current value
            Q=Q_o-y_bar*eye(size(Q_o,1));
            F=F_o;
            R=R_o-y_bar*eye(size(R_o,1));             % [Q,F;F',R]=S(lambda)-y*I
        end
        
        % ------------------------------------------- %
        %   Step 2: check if R is negative definite   %    % R<0 is a necessary condition !!
        % ------------------------------------------- %
        
        max_R_eig=IS_R_NEG(R,nth_lmi,nA);
        
        if more_constraint==1,
            if feasible,
                yh_u=y_ubound;
                yh_l=y_bar;
                dist=abs(y_ubound-y_bar);
                beta=0.5;
                while abs(yh_u-yh_l)>=0.05*dist
                    yh_m=beta*yh_l+(1-beta)*yh_u;
                    R=R_o+yh_m*R_r;
                    max_R_eig=IS_R_NEG(R,nth_lmi,nA);
                    if more_constraint==0,
                        yh_u=yh_m;
                    elseif more_constraint==1,
                        yh_l=yh_m;
                    end
                end
                
                % re-compute Q,F,R by the new testing point
                y_bar=yh_u;
                Q=Q_o+y_bar*Q_r;
                F=F_o+y_bar*F_r;
                R=R_o+y_bar*R_r;
                more_constraint=0;
            else
                y_bar=y_bar+(max_R_eig+1e-6);            % move y_bar back!
                Q=Q-(max_R_eig+1e-6)*eye(size(Q,1));
                R=R-(max_R_eig+1e-6)*eye(size(R,1));
                more_constraint=0;
            end
        end
        
        % ---------------------------------------------------------------- %
        %   Step 3: <<R is negative definite>> check if the Hamiltonian    %
        %   matrix has eigenvalues on the imaginary axis                   %
        % ---------------------------------------------------------------- %
        
        % -------------------- %
        %  Hamiltonian matrix  %
        % -------------------- %
        
        H=Hamiltonian_matrix(A,B,Q,F,R);
        IS_Hsys_OK(H,A,B,Q,F,R,nth_lmi);
        
        if more_constraint==1,
            yh_u=y_ubound;
            yh_l=y_bar;
            dist=abs(y_ubound-y_bar);
            beta=0.5;
            
            while abs(yh_u-yh_l)>=0.05*dist
                yh_m=beta*yh_l+(1-beta)*yh_u;
                if feasible,
                    Q=Q_o+yh_m*Q_r;
                    F=F_o+yh_m*F_r;
                    R=R_o+yh_m*R_r;
                else
                    Q=Q_o-yh_m*eye(size(Q,1));
                    R=R_o-yh_m*eye(size(R,1));
                end
                H=Hamiltonian_matrix(A,B,Q,F,R);
                IS_Hsys_OK(H,A,B,Q,F,R,nth_lmi);
                if more_constraint==0,
                    yh_u=yh_m;
                elseif more_constraint==1,
                    yh_l=yh_m;
                end
            end
            
            % get a new y upper bound !
            y_ubound=yh_u;
            more_constraint=0;
        end
    end   %%% <--- end of for loop of 'fdeplmi' %%%
    
    
    %% ================== %%
    %%                    %%
    %%  Stoping criteria  %%
    %%                    %%
    %% ================== %%
    
    if feasible==0,  % status: Phase 1
        if no_improve>5
            if y_ubound<0,
                feasible=1;
                re_solve_LP=1;
                L.xf=lambda;
                if ~isempty(L.obj)    % mincx problem: switch to phase 2
                    [LP_A,LP_B,estimate_obj]=prep_PS2(L,lambda,ndecv_o,phc_A,phc_B);
                    y_ubound=0;
                    counter=0;
                    no_improve=1;
                else
                    done=1;         % feasibility problem
                end
            end
            
        elseif re_solve_LP==0,
            y_ubound=y_bar;
            if y_ubound<0 && abs(y_ubound-y)<1e-2
                feasible=1;
                re_solve_LP=1;
                L.xf=lambda;
                if ~isempty(L.obj)      % mincx problem: switch to phase 2
                    [LP_A,LP_B,estimate_obj]=prep_PS2(L,lambda,ndecv_o,phc_A,phc_B);
                    y_ubound=0;
                    counter=0;
                    no_improve=1;
                else
                    done=1;             % feasibility problem
                end
                
            elseif y_ubound<0 && abs(y_ubound-y)>=1e-2 % y_ubound is not satisfactory
                no_improve=1;                          % reset no_improve counter and
                % keep working ...
            elseif y_ubound>0 && abs(y_ubound-y)<1e-2  % infeasibility established
                gain=[];
                xopt=[];
                disp(' ')
                disp('The problem is infeasible !')
                return
            end
            
        elseif re_solve_LP==1,  % there is no improve on y_ubound
            if y_ubound<0,
                no_improve=no_improve+1; % start no_improve counter
            end
        end
        
        %%% display ... %%%
        if counter~=0
            str0=num2str(y_ubound);
            str1=[blanks(2),sft(num2str(counter),5),blanks(17),str0];
            disp(str1)
        end
        
    elseif feasible==1,    % status: Phase 2
        if re_solve_LP==0,
            y_ubound=y_bar;
            estimate_obj=(L.obj*[lambda;1]+(y_ubound)*L.obj*[L.xf;1])/(1+y_ubound);
            no_improve=1;
            if abs(y_ubound-y)<error_tol2
                done=1;
            end
        elseif re_solve_LP==1 && abs(y_ubound-y)<error_tol1,
            if no_improve>=5,
                done=1;
            else
                no_improve=no_improve+1;
            end
        end
        
        %%% display ... %%%
        str2=[blanks(2),sft(num2str(counter),5),blanks(17),num2str(estimate_obj)];
        disp(str2)
    end
    
    counter=counter+1;  % increase counter
    
end %%% <--- end of Mail Loop!!

if feasible==1
    if ~isempty(L.obj)    % mincx problem !
        xopt=(lambda+(y_ubound)*L.xf)/(1+y_ubound);
        cost=L.obj*[xopt;1];
    else                  % feasibility problem !
        xopt=L.xf;
        cost=[];
        disp(' ')
        disp('Problem is feasible ... A feasible point is found!')
    end
else
    cost=[];
    xopt=[];
    disp(' ')
    disp(['No feasible point found after ',num2str(counter-1),' iteration'])
end

global ABST
ABST.E=E;
ABST.xopt=xopt;
fdlmi_value;


%%%%%%%%%%%% Inner Functions %%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
%%                  %%
%%  Initialization  %%
%%                  %%
%%%%%%%%%%%%%%%%%%%%%%

function y_ubound=initialize(nfdeplmi,fdeplmi)

global L LP_A LP_B
global ndecv_e ndecv_o %#ok<*NUSED>
global more_constraint
global feasible

lmi_sz_m1=[nfdeplmi;fdeplmi];
nlmi=length(lmi_sz_m1);
count=0;
for flcnt=1:ndecv_o
    nth_lmi=lmi_sz_m1(count+1);
    sz = L.lmi{nth_lmi}.size;
    A  = L.lmi{nth_lmi}.A;
    nA = size(A,1);
    vp = [zeros(nA,1);rand(sz-nA,1)];
    [coeff,constant]=hyperplane(L,nth_lmi,vp,ndecv_o);
    if feasible,
        coeff_y=vp'*L.lmi{nth_lmi}.MF*vp;
    else
        coeff_y=-(vp'*vp);
    end
    cnst=[coeff,coeff_y];
    norm_cnst=norm(cnst,2);
    cnst=cnst/norm_cnst;    % normalize the constraint
    LP_A=[LP_A;cnst];       % incorperate the constraint with previous constraints
    LP_B=[LP_B;-constant/norm_cnst];
    count=count+1;
    if (count+1)>nlmi,
        count=0;
    end
end

eigc=[];
for flcnt=1:nlmi
    nth_lmi=lmi_sz_m1(flcnt);
    CST = L.lmi{nth_lmi}.C;
    sz  = L.lmi{nth_lmi}.size;
    A   = L.lmi{nth_lmi}.A;
    nA  = size(A,1);
    CST = CST(nA+1:sz,nA+1:sz);
    eigc=[eigc;eig(0.5*(CST+CST'))];
end
y_ubound=max(eigc);


%%%%%%%%%%%%%%%%
%%            %%
%%  Solve LP  %%
%%            %%
%%%%%%%%%%%%%%%%

function [lambda,y]=LP_solve(LP_C,method)

global L LP_A LP_B

switch method
    case 'a'
        [lambda_ext,obj,status,popt]=lps(-LP_A,-LP_B,LP_C,[1e-6,1e-12,1e-12,1000]); %#ok<*ASGLU>
        if status==2,
            error('*** The problem is not feasible!! ***')
        elseif status==3,
            LP_A=[[zeros(1,ndecv_o),-1];LP_A];
            LP_B=[1;LP_B];
            [lambda_ext,obj,status,popt]=lps(-LP_A,-LP_B,LP_C,[1e-6,1e-12,1e-12,1000]);
        end
    case 'm',
        [lambda_ext,multiplier,how]=lp(LP_C,LP_A,LP_B);
        if strcmp(how,'infeasible')
            error('*** IQC is not feasible!! ***')
        elseif strcmp(how,'unbounded')
            LP_A=[[zeros(1,ndecv_o),-1];LP_A];
            LP_B=[1;LP_B];
            [lambda_ext,multiplier,how]=lp(LP_C,LP_A,LP_B);
        end
end

lambda  = lambda_ext(1:length(lambda_ext)-1);
y       = lambda_ext(length(lambda_ext));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              %%
%%  Verify NFDLMIs and compute  %%
%%  the polyhedral constraints  %%
%%                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function max_MX_eig=V_NFDLMI(MX,nth_lmi)

global L LP_A LP_B
global feasible
global ndecv_o
global more_constraint
global re_solve_LP

max_MX_eig=[];
[cholR,cholp]=chol(-MX);
if cholp>0                                  % there are positive eigenvalues !
    % ------------------------- %
    %   computing constraints   %
    % ------------------------- %
    [V,D]=eig(MX);
    pos=find(diag(D)>-1e-6);
    for flcnt2=1:length(pos)
        vp=V(:,pos(flcnt2));
        [coeff,constant]=hyperplane(L,nth_lmi,vp,ndecv_o);
        if feasible,
            coeff_y=vp'*L.lmi{nth_lmi}.MF*vp;
        else
            coeff_y=-(vp'*vp);
        end
        cnst=[coeff,coeff_y];
        norm_cnst=norm(cnst,2);
        cnst=cnst/norm_cnst;       % normalize the constraint
        LP_A=[LP_A;cnst];          % incorperate the constraint with previous constraints
        LP_B=[LP_B;-constant/norm_cnst];
    end
    max_MX_eig=max(diag(D));
    more_constraint=1;         % set flag to 1
    re_solve_LP=1;
else
    more_constraint=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         %%
%%  Is R sign definite ?   %%
%%                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function max_R_eig=IS_R_NEG(R,nth_lmi,nA)

global L LP_A LP_B
global ndecv_o
global feasible
global more_constraint
global re_slove_LP

max_R_eig=[];

[cholR,cholp]=chol(-R);
if cholp>0                                   % there are positive eigenvalues !
    % ------------------------- %
    %   computing constraints   %
    % ------------------------- %
    [V,D]=eig(R);
    pos=find(diag(D)>-1e-6);
    for flcnt2=1:length(pos)
        vp=V(:,pos(flcnt2));
        vp_r=[zeros(nA,1);vp];
        [coeff,constant]=hyperplane(L,nth_lmi,vp_r,ndecv_o);
        if feasible,
            coeff_y=vp_r'*L.lmi{nth_lmi}.MF*vp_r;
        else
            coeff_y=-(vp_r'*vp_r);
        end
        cnst=[coeff,coeff_y];
        norm_cnst=norm(cnst,2);
        cnst=cnst/norm_cnst;       % normalize the constraint
        LP_A=[LP_A;cnst];          % incorperate the constraint with previous constraints
        LP_B=[LP_B;-constant/norm_cnst];
    end
    more_constraint=1;         % set flag to 1
    max_R_eig=max(diag(D));
    re_solve_LP=1;
else
    more_constraint=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              %%
%%  Compute Hamiltonian matrix  %%
%%                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H=Hamiltonian_matrix(A,B,Q,F,R)

H11=A-B*inv(R)*F'; %#ok<*MINV>
H12=B*inv(R)*B';
H21=Q-F*inv(R)*F';
H22=-A'+F*inv(R)*B';
H=[H11,H12;H21,H22];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                      %%
%%  Hamiltonian system has eigenvalues  %%
%%  on the imaginary axis               %%
%%                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IS_Hsys_OK(H,A,B,Q,F,R,nth_lmi)

global L LP_A LP_B
global ndecv_o
global feasible
global more_constraint
global re_solve_LP

[VH,DH]=eig(H);
pos=find(diag(abs(real(DH)))<1e-6);
imag_DH=diag(imag(DH));     % imaginery part of eigenvalues
imag_VH=VH(:,pos);          % eig-vectors corresponding eigenvalues
% which are on the imaginery axis
freq_V=imag_VH(:,imag_DH(pos)>=0);

if ~isempty(pos)
    diag(DH);
    omega=imag_DH(pos);
    pos= imag_DH(pos)>=0;
    omega=omega(pos);
    
    for flcnt10=1:length(omega)
        oga=omega(flcnt10);
        GG1=[inv(1i*oga*eye(size(A,1))-A)*B;eye(size(B,2))];
        GG=GG1'*[Q,F;F',R]*GG1;
        [VGG,DGG]=eig(GG);
        eig_r=real(diag(DGG));
        pos=(find(eig_r>-1e-6));
        % eig_r(pos)
        for flcnt2=1:length(pos)
            VGG_pos=VGG(:,pos(flcnt2));
            vec=GG1*VGG_pos;
            [coeff,constant]=hyperplane(L,nth_lmi,vec,ndecv_o);
            coeff=real(coeff);
            constant=real(constant);
            if feasible,
                coeff_y=real(vec'*L.lmi{nth_lmi}.MF*vec);
            else
                coeff_y=-real(vec'*vec);
            end
            cnst=[coeff,coeff_y];
            norm_cnst=norm(cnst,2);
            cnst=cnst/norm_cnst;       % normalize the constraint
            LP_A=[LP_A;cnst];          % incorperate the constraint with previous constraints
            LP_B=[LP_B;-constant/norm_cnst];
        end
    end
    more_constraint=1;         % set flag to 1
    re_solve_LP=1;
else
    more_constraint=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            %%
%%  prepare necessary change  %%
%%  for Phase 2               %%
%%                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LP_A,LP_B,estimate_obj]=prep_PS2(L,lambda,ndecv_o,phc_A,phc_B)
% function [L,LP_A,LP_B,estimate_obj]=prep_PS2(L,lambda,ndecv_o,phc_A,phc_B)
% Prepare the necessary change to Phase 2
nobj=length(L.obj);
LP_A=[zeros(1,ndecv_o),-1;L.obj(1:nobj-1),-1];
LP_B=[1;L.obj(1:nobj-1)*lambda+L.obj(nobj)];

if ~isempty(phc_A)
    AA=phc_A(:,1:ndecv_o);
    phc_A=[AA,(AA*lambda-phc_B)];
    LP_A=[LP_A;phc_A];
    LP_B=[LP_B;phc_B];
end

estimate_obj=L.obj(1:nobj-1)*lambda+L.obj(nobj);
disp(' ')
disp('***** feasible point found! Switch to P2 ...')
disp(' ')
%%             1         2         3         4
%% ## 1234567890123456789012345678901234567890
disp(' Iteration    :      estimate objective')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   %%
%%  compute ployhegral constraints   %%
%%                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coeff,constant]=hyperplane(L,nth_lmi,vp,ndecv_o)
% compute polyhedral constraint
nterm=size(L.lmi{nth_lmi}.X,2);
coeff=zeros(1,ndecv_o);              % compute constraint
constant=0;
for flcnt3=1:nterm
    G=L.lmi{nth_lmi}.G{flcnt3}';
    X=L.lmi{nth_lmi}.X{flcnt3};
    H=L.lmi{nth_lmi}.H{flcnt3};
    if ~isempty(X)
        coeff=coeff+form_constraint(G,H,X,vp,ndecv_o);      % non-constant term
    else
        constant=constant+vp'*L.lmi{nth_lmi}.C*vp;  % constant term
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   %%
%%  compute the coefficients of      %%
%%  polyhedral constraints           %%
%%                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coeff]=form_constraint(A,B,X,v,nvar)
% function [coeff]=form_constraint(A,B,X,v,nvar)
%
% given v, and matrix function A'*X*B, where X is a matrix
% variable; this program computes the coefficients of the
% hyperplane H(x)=(v'A')*X*(Bv)
%
% written by cykao@mit.edu  Oct. 14 1998
% last modified Apr. 6 1999

mX=size(X,1);
nX=size(X,2);
v_left=A*v;
v_right=B*v;

coeff=zeros(1,nvar);
if (mX==1 && nX==1)
    nth_var=X(1,1);
    if nth_var>0
        coeff(nth_var)=v_left'*v_right;
    elseif nth_var<0
        coeff(-nth_var)=-v_left'*v_right;
    end
else
    for flcnt1=1:mX                    % index of row of matrix X
        for flcnt2=1:nX                % index of column of matrix X
            nth_var=X(flcnt1,flcnt2);  % index of independent variables
            if nth_var~=0,
                if nth_var>nvar,
                    error('Bad index of independent variables!!')
                else
                    if nth_var>0
                        coeff(nth_var)=coeff(nth_var)+conj(v_left(flcnt1))*v_right(flcnt2);
                    elseif nth_var<0
                        coeff(-nth_var)=coeff(-nth_var)-conj(v_left(flcnt1))*v_right(flcnt2);
                    end
                end
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%                    %%
%%  plug in \lambda   %%
%%                    %%
%%%%%%%%%%%%%%%%%%%%%%%%

function MX=lambda_eval(lambda,nth_lmi)

global L

nterm=length(L.lmi{nth_lmi}.X);
MX=zeros(L.lmi{nth_lmi}.size);
for flcnt2=1:nterm
    G = L.lmi{nth_lmi}.G{flcnt2};
    X = L.lmi{nth_lmi}.X{flcnt2};
    H = L.lmi{nth_lmi}.H{flcnt2};
    if ~isempty(X)
        MX= MX + matrix_eval(G,X,H,lambda);  % non-constant term in LMI
    else
        MX= MX + L.lmi{nth_lmi}.C;           % constant term in LMI
    end
end
MX=0.5*(MX+MX');   %%% all LMIs "M<0" are interpreted as " M+M'<0 " %%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                 %%
%%  matrix variable evaluation     %%
%%                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xe]=matrix_eval(A,X,B,coeff)
% function [Xe]=matrix_eval(A,X,B,coeff)
%
% evaluate matric function F(x)=A*X*B
%
% written by cykao@mit.edu  Oct. 14 1998
% last modified  Apr. 6 1999

mX=size(X,1);
nX=size(X,2);
lc=length(coeff);

Xv=zeros(mX,nX);

for flcnt1=1:mX      % index of row of matrix X
    for flcnt2=1:nX  % index of column of matrix X
        nth_var=X(flcnt1,flcnt2);
        if nth_var>lc,
            error('Bad index of independent variables')
        else
            if nth_var>0,
                Xv(flcnt1,flcnt2)=coeff(nth_var);
            elseif nth_var<0
                Xv(flcnt1,flcnt2)=-coeff(-nth_var);
            end
        end
    end
end
Xe=A*Xv*B;
