function gain=iqc_analysis_l2gain_lmilab(f,z)
% function gain=iqc_analysis_l2gain_lmilab(f,z)
%
% gives best estimate of the L2 gain f-> z in the system of IQC's
% defined using the "abst" environment "iqc"
%
% stores the optimal multipliers in global variable ABSTSOLUTION
% (so that they can be retrieved using "value")
%
% 'f' and 'z' can be of type "inp" or "sgn"
%
% The default SDP solver is the solver coming with the LMI Control Toolbox
% by The MathWorks, Inc.
%
% if ABST.lmiparameter(5)>100, the resulting LMI script is
% exported in the format of two files: *_exe.m and *_exe.mat
% (containing the script and the datafile respectively)
%
% Written by ameg@mit.edu, last modified November 10, 1997
% Modified by cykao@mit.edu on Oct. 25, 1998
% Last modified by cykao@ee.mu.oz.au on Aug. 14 2004
% Last modified by cmj on 2013/4/23

global ABST
str={};

ABST.systemtype = 'continuous';

lmiparameter=ABST.lmiparameter;

%vrb=0 to suppress the output
vrb=(lmiparameter(5)==0)||(lmiparameter(5)==777);

% first, process the iqc abst log
E=iqc_extract(f,z);

% This line was added by C Kao on Jun. 15 1998 for 'link' stuff
% --------->
E=iqc_reduce(E);

% symbolic names for interior types:
var=2;
lmi=4;


if vrb,
    if lmiparameter(5)>100,
        disp_str(50,'iqc_analysis_l2gain_lmilab')
    else
        disp_str(51,'iqc_analysis_l2gain_lmilab')
    end
end

export_date=date;
str{1}=['%% Created by iqc_analysis_l2gain_lmilab on ',export_date];
str{2}='\n%% Load coefficient matrices ...';
str{3}='load iqc_analysis_l2gain_lmilab_exe';

%% initialize the LMI Control Toolbox
str{4}='\n%% Initialize the LMI Lab';
str{5}='setlmis([]);';
eval(str{5});

%% variable definite
if vrb,
    disp_str(52)
end
str{6}='\n%% Define multiplier variables ...';
sc=7; % string counter

kvar=find(E.T==var);
for k=1:length(kvar);
    ks=num2str(k);
    eval(['struct' ks '=E.L{' num2str(kvar(k)) '};'])
    str{sc}=['x' ks '=lmivar(3,struct' ks ');'];
    eval(str{sc});
    sc=sc+1;
end

%% now we define the "non-KYP" lmi
if vrb,
    disp_str(53)
end
str{sc}='\n%% Define non-KYP lmis ...';
sc=sc+1;
kvar=find(E.T==lmi);
for k=1:length(kvar),
    ks=num2str(k);
    kk=kvar(k);
    [nn,mm]=size(E.AB{kk});    % dimensions
    nns=num2str(nn);
    vs=mm-nn;                       % LMI size
    vss=num2str(vs);
    if nn>0,
        eval(['state' ks '=E.AB{kk};']);
        str{sc}=['p',ks,'=lmivar(1,[',nns,' 1]);']; %#ok<*AGROW>
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
    for i=1:size(E.X{kk},2),
        is=num2str(i);
        nvs=num2str(E.X{kk}(1,i));   % variable number
        eval(['l_' ks '_' is '=E.C{kk}(ng+1:ng+E.X{kk}(2,i),:)'';'])
        eval(['r_' ks '_' is '=E.C{kk}(nh+1:nh+E.X{kk}(3,i),:);'])
        str{sc}=['lmiterm([',ks,' 1 1 ',nvs,'],l_',ks,'_',is,...
            ',r_',ks,'_',is,',''s'');']; %#ok<AGROW>
        eval(str{sc});
        sc=sc+1;
        nh=nh+E.X{kk}(3,i);
        ng=ng+E.X{kk}(2,i);
    end
end

%% now we define the "main" lmi
if vrb,
    disp_str(54)
end
str{sc}='\n%% KYP LMI terms...';
sc=sc+1;
k=E.nlmi+1;
ks=num2str(k);
nn=E.nstates;
nns=num2str(nn);
vs=E.ninputs;
vss=num2str(vs);
mm=nn+vs; %#ok<*NASGU>

if nn>0,
    str{sc}=['[P,nP,sP]=lmivar(1,[',nns,' 1]);'];
    eval(str{sc});
    sc=sc+1;
    main_state=E.ab;
    if strcmp(ABST.systemtype,'continuous'),
        str{sc}=['lmiterm([',ks,' 1 1 P],[eye(',...
            nns,');zeros(',vss,',',nns,')],main_state,''s'');'];
        eval(str{sc});
        sc=sc+1;
    elseif strcmp(ABST.systemtype,'discrete'),
        str{sc}=['lmiterm([',ks,' 1 1 P],main_state'',main_state);'];
        eval(str{sc});
        sc=sc+1;
        str{sc}=['lmiterm([',ks,' 1 1 P],[-eye(',nns,') zeros(',nns,',',vss,...
            ')]'',[eye(',nns,') zeros(',nns,',',vss,')]);'];
        eval(str{sc});
        sc=sc+1;
    else
        disp_str(55)
    end
end
for j=1:E.nsimple
    js=num2str(j);
    nh=0;
    ng=0;
    n=E.simples(j);
    for i=1:size(E.X{n},2)
        is=num2str(i);
        vrs=num2str(E.X{n}(1,i));
        eval(['L_' js '_' is ...
            '=E.B{E.simples(j)}(ng+1:ng+E.X{E.simples(j)}(2,i),:)'';'])
        eval(['R_' js '_' is ...
            '=E.gqfm(j)*E.C{E.simples(j)}(nh+1:nh+' ...
            'E.X{E.simples(j)}(3,i),:);'])
        str{sc}=['lmiterm([',ks,' 1 1 ',vrs,'],L_',js,'_',is,...
            ',R_',js,'_',is,',''s'');'];
        eval(str{sc});
        sc=sc+1;
        nh=nh+E.X{E.simples(j)}(3,i);
        ng=ng+E.X{E.simples(j)}(2,i);
    end
end

% introduce the gain estimation terms
% -----------------------------------
str{sc}='\n%% L2 gain estimation terms ...';
sc=sc+1;
c_z=E.C{double(z)};
c_f=E.C{double(f)};
str{sc}=['lmiterm([',ks,' 1 1 0],c_z''*c_z);'];
eval(str{sc});
sc=sc+1;
str{sc}='minimize_variable=lmivar(1,[1 1]);';
eval(str{sc});
sc=sc+1;
str{sc}=['lmiterm([-',ks,' 1 1 minimize_variable],c_f''*c_f,1);'];
eval(str{sc});
sc=sc+1;
k=k+1;
ks=num2str(k);

%% solving the system of LMI's
str{sc}='\n%% Solving the system of LMIs ...';
sc=sc+1;
str{sc}='lmi=getlmis;';
eval(str{sc});
sc=sc+1;
str{sc}='ndec=decnbr(lmi);';
eval(str{sc});
sc=sc+1;

if vrb
    disp_str(56,num2str(ndec))
end

str{sc}='c=zeros(ndec,1);';
eval(str{sc});
sc=sc+1;
str{sc}='nc=decinfo(lmi,minimize_variable);';
eval(str{sc});
sc=sc+1;
str{sc}='c(nc)=1;';
eval(str{sc});
sc=sc+1;
gain = 0;
str{sc}='[gain,xopt]=mincx(lmi,c,lmiparameter);';
eval(str{sc});
sc=sc+1;
if ~isempty(gain) && gain > 0
    if vrb
        disp_str(57,'L2-gain','Feasible')
    end
    str{sc}='gain=sqrt(gain)';
    eval([str{sc},';']);
    sc=sc+1;
else
    if vrb
        disp_str(57,'L2-gain','Infeasible')
    end
    str{sc}='gain=inf';
    eval([str{sc},';']);
    sc=sc+1;
end


%%
if lmiparameter(5)>100;
    save iqc_analysis_l2gain_lmilab_exe main_state lmiparameter c_f c_z
    if exist('state1','var')
        save iqc_analysis_l2gain_lmilab_exe state* -append
    end
    if exist('struct1','var');
        save iqc_analysis_l2gain_lmilab_exe struct* -append
    end
    if exist('l_1_1','var');
        save iqc_analysis_l2gain_lmilab_exe l_* -append
    end
    if exist('r_1_1','var');
        save iqc_analysis_l2gain_lmilab_exe r_* -append
    end
    if exist('L_1_1','var');
        save iqc_analysis_l2gain_lmilab_exe L_* -append
    end
    if exist('R_1_1','var');
        save iqc_analysis_l2gain_lmilab_exe R_* -append
    end
    fid=fopen('iqc_analysis_l2gain_lmilab_exe.m','wt');
    for i=1:sc-1
        fprintf(fid,[str{i} '\n']);
    end
    fclose(fid);
end

ABST.xopt=xopt;
ABST.E=E;
ABST.c_f=c_f;
ABST.c_z=c_z;

if ~isempty(xopt);
    ABST.P=decs2mats(xopt,sP);
    ABST.xg=decs2mats(xopt,nc);
    iqc_value;
end

%%
function [E]=iqc_reduce(E)
%% function [E]=iqc_reduce(E)
%
% -------------------------------------------------
% the following function is for reducing the system
% program done by C. Kao on Jun. 15 1998
% last modified on June 16 1998
% last modified on June 17 1998
% last modified on June 21 1998
% --------------------------------------------------
% internal function:
% the purpose of this internal function is dealing with
% the 'link' stuff and reducing the system

global ABST
A=ABST;

% symbolic names for interior types:
% -----------------------------------
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

% collect the log index in which 'C' field has been used
% -------------------------------------------------------
lnk_index=find(E.T==lnk);
inp_index=find(E.T==inp);
sgn_index=find(E.T==sgn);
vsg_index=find(E.T==vsg);
csg_index=find(E.T==csg);
vcs_index=find(E.T==vcs);
qfm_index=find(E.T==qfm);
c_index=[inp_index;sgn_index;vsg_index;csg_index;vcs_index;qfm_index];

% collect the log index in which the 'B' field contain the
% output coeffiecients.
% ---------------------------------------------------------
b_index=qfm_index;

% define some useful constants
% ------------------------------
ns=E.nstates;        % number of total states
ni=E.ninputs;        % number of total inputs
AB=E.ab;             % state space matrix [A,B]
[mAB,nAB]=size(AB);  % size of [A,B]
BigC=[];             % BigC will be used to collect all non-empty entries in C field
BigB=[];             % BigB will be used to collect the entries of B field which
% correspond to output coefficients.

% collect all inputs which are linked to signals
% Form 'Cx+Dw=0'
% ----------------------------------------------
lnk_inputs=[];
CD_lnk=[];
for forlp_counter=1:length(lnk_index)
    a1=A.log(lnk_index(forlp_counter),6);
    a2=A.log(lnk_index(forlp_counter),7);
    slnk=num2str(lnk_index(forlp_counter));
    err_msg1=['Opps ... Error about ''link'' occured in log #',slnk];
    
    if A.log(a1,1)==inp && A.log(a2,1)==inp,
        error('Opps! something wrong in iqc_extract.m; input==input found!')
    elseif A.log(a1,1)==inp,
        lnk_inputs=[lnk_inputs;a1];
    elseif A.log(a2,1)==inp,
        lnk_inputs=[lnk_inputs;a2];
    else
        error([err_msg1,' ,''link'' without inputs'])
    end
    
    CD_lnk=[CD_lnk;E.C{lnk_index(forlp_counter)}];
end
[mCD_lnk,nCD_lnk]=size(CD_lnk);

% check if any input linked to more than one signal
% --------------------------------------------------
for forlp_counter1=1:length(lnk_index)-1
    for forlp_counter2=forlp_counter1+1:length(lnk_index)
        if lnk_inputs(forlp_counter1)==lnk_inputs(forlp_counter2)
            error('There is an input being linked to more than one signals')
        end
    end
end

% collecting all output coeffiecents in C field and B field together
% ------------------------------------------------------------------
c_index_new=[];
for forlp_counter=1:length(c_index)
    if  c_index(forlp_counter)<=size(E.C,2)
        if ~isempty(E.C{c_index(forlp_counter)})
            %             c_index(forlp_counter);
            BigC=[BigC;E.C{c_index(forlp_counter)}];
            c_index_new=[c_index_new;c_index(forlp_counter)];
        end
    end
end
[mBigC,nBigC]=size(BigC); %#ok<*ASGLU>

b_index_new=[];
for forlp_counter=1:length(b_index)
    if  b_index(forlp_counter)<=size(E.B,2)
        if ~isempty(E.B{b_index(forlp_counter)})
            b_index(forlp_counter);
            BigB=[BigB;E.B{b_index(forlp_counter)}];
            b_index_new=[b_index_new;b_index(forlp_counter)];
        end
    end
end
[mBigB,nBigB]=size(BigB);

% concatenate BigB and BigC matrices
BigBC=[BigB;BigC];
[mBigBC,nBigBC]=size(BigBC);



% main loop: recomputing the state space matrix [A,B] and the output
%            coefficient matrix
% -------------------------------------------------------------------

% safety check
% ------------
if length(lnk_inputs)~=length(lnk_index)
    error('An unexpected error! ... in iqc_reduce.m')
end

elim=[];                    % elim will be used to collect the inputs which
% will be removed (the linked inputs)
nreplaced=0;
cnt=0;
for forlp_counter=1:length(lnk_inputs)
    inp_replaced=lnk_inputs(forlp_counter);
    vs=A.log(lnk_index(forlp_counter),2);
    hs=A.log(lnk_index(forlp_counter),3);
    nreplaced=nreplaced+vs;                 % counting how many signals and
    % inputs are to be replaced
    
    pos1=E.POS(inp_replaced):E.POS(inp_replaced)+vs-1;
    pos2=[ns+1:E.POS(inp_replaced)-1,E.POS(inp_replaced)+vs:ns+ni];
    elim=[elim,pos1];
    
    % check if this link has been processed before
    chk1=CD_lnk(cnt+1:cnt+vs,:);
    chk2=(chk1==zeros(size(chk1,1),size(chk1,2)));
    if ~all(chk2(:))
        D_replaced_lnk=CD_lnk(cnt+1:cnt+vs,pos1);
        r=rank(D_replaced_lnk,1e-6);
        if r~=vs,
            error([err_msg1,' ,non-invertable matrix found!'])
        else
            invD_replaced_lnk=inv(D_replaced_lnk);
            D_lnk=invD_replaced_lnk*(-1)*CD_lnk(cnt+1:cnt+vs,pos2); %#ok<*MINV>
            C_lnk=invD_replaced_lnk*(-1)*CD_lnk(cnt+1:cnt+vs,1:ns);
        end
        
        % re-compute state space matrix
        % -----------------------------
        B_inp1=AB(:,pos1);
        B_reduced=AB(:,pos2)+B_inp1*D_lnk;
        A_reduced=AB(:,1:ns)+B_inp1*C_lnk;
        B_reduced_ext=zeros(mAB,nAB-ns);
        B_reduced_ext(:,pos2-ns)=B_reduced;
        AB=[A_reduced,B_reduced_ext];
        
        % re-compute output coefficient matrix
        % ------------------------------------
        D_inp1=BigBC(:,pos1);
        D_reduced=BigBC(:,pos2)+D_inp1*D_lnk;
        C_reduced=BigBC(:,1:ns)+D_inp1*C_lnk;
        D_reduced_ext=zeros(mBigBC,nBigBC-ns);
        D_reduced_ext(:,pos2-ns)=D_reduced;
        BigBC=[C_reduced,D_reduced_ext];
        
        % re-compute the [C,D] matrix of 'Cx+Dw=0'
        % ----------------------------------------
        D_inp1_CD_lnk=CD_lnk(:,pos1);
        D_reduced_CD_lnk=CD_lnk(:,pos2)+D_inp1_CD_lnk*D_lnk;
        C_reduced_CD_lnk=CD_lnk(:,1:ns)+D_inp1_CD_lnk*C_lnk;
        D_reduced_ext_CD_lnk=zeros(mCD_lnk,nCD_lnk-ns);
        D_reduced_ext_CD_lnk(:,pos2-ns)=D_reduced_CD_lnk;
        CD_lnk=[C_reduced_CD_lnk,D_reduced_ext_CD_lnk];
    end
    
    % increase counter
    cnt=cnt+vs;
end


% build the index matrix which will be used to re-construct the
% state space matrix [A,B] and output matrices in C field
% --------------------------------------------------------------
elim=sort(elim);
if isempty(elim)
    re_construct=1:nBigC;
else
    re_construct=[];
    cnt=0;
    for forlp_counter=1:length(elim)
        if forlp_counter==length(elim)
            re_construct=[re_construct,cnt+1:elim(forlp_counter)-1];
            re_construct=[re_construct,elim(forlp_counter)+1:nBigC];
        else
            re_construct=[re_construct,cnt+1:elim(forlp_counter)-1];
            cnt=elim(forlp_counter);
        end
    end
end

% re-construct the state space matrix [A,B] and output matrices
% in C field
% ---------------------------------------------------------------
E.ab=AB(:,re_construct);
BigB=BigBC(1:mBigB,re_construct);
BigC=BigBC(mBigB+1:mBigBC,re_construct);

cnt=0;
for forlp_counter=1:length(c_index_new)
    vs=size(E.C{c_index_new(forlp_counter)},1);
    E.C{c_index_new(forlp_counter)}=BigC(cnt+1:cnt+vs,:);
    cnt=cnt+vs;
end

cnt=0;
for forlp_counter=1:length(b_index_new)
    vs=size(E.B{b_index_new(forlp_counter)},1);
    E.B{b_index_new(forlp_counter)}=BigB(cnt+1:cnt+vs,:);
    cnt=cnt+vs;
end

cnt=0;
for forlp_counter=1:length(lnk_index)
    vs=size(E.C{lnk_index(forlp_counter)},1);
    E.C{lnk_index(forlp_counter)}=CD_lnk(cnt+1:cnt+vs,:);
    cnt=cnt+vs;
end

E.ninputs=E.ninputs-nreplaced;