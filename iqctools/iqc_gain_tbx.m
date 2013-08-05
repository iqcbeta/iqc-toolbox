function [gain, varargout]=iqc_gain_tbx_YALMIP(varargin)
% function gain=lmi_gain_tbx(f,z,solver_option)
%
% gives best estimate of the L2 gain f->z in the system of IQC's
% defined using the "abst" environment "iqc"
%
% stores the optimal multipliers in global variable ABSTSOLUTION
% (so that they can be retrieved using "value")
%
% 'f' and 'y' can be of type "inp" or "sgn" 
%
% 'solver_option' is an optional input which specifies which SDP solvers
% is to be used. IQC beta use YALMIP as the gateway to various SDP
% solvers. See YALMIP manual (http://control.ee.ethz.ch/~joloef/yalmip.msql) 
% for more details
%     
%
% The default SDP solver is the solver coming with the LMI Control Toolbox 
% by The MathWorks, Inc.
%
% if ABST.options(5)>100, the resulting LMI script is
% exported in the format of two files: lmi_exe.m and lmi_exe.mat
% (containing the script and the datafile respectively)
%
% Written by ameg@mit.edu, last modified November 10, 1997
% Modified by cykao@mit.edu on Oct. 25, 1998
%
% Last modified by cykao@ee.mu.oz.au on Aug. 14 2004 

%
% corrected a couple of parentheses typos.
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu


if nargin == 3,
   if ischar(varargin{3})~=1, 
      error('unknown solver specified') 
   end
   f = varargin{1};
   z = varargin{2};
   solver_option = varargin{3};
elseif nargin<3,
   if nargin<2, 
      error('two inputs required') 
   end
   f = varargin{1}; 
   z = varargin{2};
   solver_option = [];
end


% first, process the iqc abst log
E=iqc_extract(f,z);

% This line was added by C Kao on Jun. 15 1998 for 'link' stuff
% --------->
E=iqc_reduce(E);

% vrb=0 to suppress the output
vrb=(E.options(5)==0)|(E.options(5)==777);

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


if (~isa(f,'abst'))|(~isa(z,'abst')),
   error('all inputs must be from the "abst" class')
end
if (E.T(double(f))~=sgn)&(E.T(double(f))~=inp),
   error('First argument must be input or signal')
end
if (E.T(double(z))~=sgn)&(E.T(double(z))~=inp),
   error('Second argument must be input or signal')
end
if vrb,
   if E.options(5)>100,
      disp('iqc_gain_tbx  (exporting the LMI script to lmi_exe.m) ...')
   else
      disp('iqc_gain_tbx ...')
   end
end


if isempty(solver_option) %% <-------- use LMI Control Toolbox ------

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
   for k=1:length(kvar),
       ks=num2str(k);
       eval(['struct' ks '=E.L{' num2str(kvar(k)) '};'])
       str{sc}=['x' ks '=lmivar(3,struct' ks ');'];
       eval(str{sc});
       sc=sc+1;
    end
      
    % now define all LMI's
    if vrb,
       disp('  defining the non-KYP LMIs ...')
    end
    kvar=find(E.T==lmi);
    str{sc}='\n%% Define non-KYP lmis ...';
    sc=sc+1;
    for k=1:length(kvar),
        kk=kvar(k);
        [nn,mm]=size(E.AB{kk});    % dimensions
        ks=num2str(k);
        nns=num2str(nn);
        vs=mm-nn;                       % LMI size
        vss=num2str(vs);
        if nn>0,
           eval(['state' ks '=E.AB{kk};']);
           str{sc}=['p' ks '=lmivar(1,[' nns ' 1]);'];
           eval(str{sc});
           sc=sc+1;
           str{sc}=['lmiterm([' ks ' 1 1 p' ks ...
                    '],[eye(' nns ');zeros(' vss ',' nns ')],state' ks ...
                    ',''s'');'];
           eval(str{sc});
           sc=sc+1;
        end
        nh=0;
        ng=sum(E.X{kk}(3,:));
        for i=1:size(E.X{kk},2),
            is=num2str(i);
            nvs=num2str(E.X{kk}(1,i));   % variable number
            eval(['l_' ks '_' is '=E.C{kk}(ng+1:ng+E.X{kk}(2,i),:)'';'])
            eval(['r_' ks '_' is '=E.C{kk}(nh+1:nh+E.X{kk}(3,i),:);'])
            str{sc}=['lmiterm([' ks ' 1 1 ' nvs '],l_' ks ...
                     '_' is ',r_' ks '_' is ',''s'');'];
            eval(str{sc});
            sc=sc+1;
            nh=nh+E.X{kk}(3,i);
            ng=ng+E.X{kk}(2,i);
        end
    end

    % now we define the "main" lmi 
    if vrb,
       disp('  defining the KYP LMIs ...')
    end
    str{sc}='\n%% KYP LMI terms...';
    sc=sc+1;
    nm=E.nlmi+1;
    nms=num2str(nm);
    mm=E.ninputs;
    mms=num2str(mm);
    if E.nstates>0,
       nss=num2str(E.nstates);
       str{sc}=['[P,ndec_p,XP]=lmivar(1,[' nss ' 1]);'];
       eval(str{sc});
       sc=sc+1;
       state0=E.ab;
       % lmiterm([-nm-1 1 1 P],1,1);  % ?????????
       str{sc}=['lmiterm([' nms ' 1 1 P],[eye(' nss ...
                ');zeros(' mms ',' nss ')],state0,''s'');'];
       eval(str{sc});
       sc=sc+1;      
    end
    for j=1:E.nsimple,
        js=num2str(j);
        nh=0;           % C-position counter for the right (C) terms 
        ng=0;           % C-position counter for the left (B) terms
        for i=1:size(E.X{E.simples(j)},2),
            is=num2str(i);
            vrs=num2str(E.X{E.simples(j)}(1,i));
            eval(['L_' js '_' is ...
                  '=E.B{E.simples(j)}(ng+1:ng+E.X{E.simples(j)}(2,i),:)'';'])
            eval(['R_' js '_' is ...
                  '=E.gqfm(j)*E.C{E.simples(j)}(nh+1:nh+' ...
                  'E.X{E.simples(j)}(3,i),:);']) 
            str{sc}=['lmiterm([' nms ' 1 1 ' vrs '],L_' js '_' is ...
                     ',R_' js '_' is ',''s'');'];
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
    str{sc}=['lmiterm([' nms ' 1 1 0],c_z''*c_z)'];
    eval(str{sc});
    sc=sc+1;
    str{sc}='xg=lmivar(1,[1 1]);';
    eval(str{sc});
    sc=sc+1;
    str{sc}=['lmiterm([-' nms ' 1 1 xg],c_f''*c_f,1);'];
    eval(str{sc});
    sc=sc+1;
 
    % solving the system of LMI's 
    str{sc}='\n%% Solving the system of LMIs ...';
    sc=sc+1;
    str{sc}='lmi=getlmis;';
    eval(str{sc});
    sc=sc+1;
    str{sc}='ndec=decnbr(lmi);';
    eval(str{sc});
    sc=sc+1;
    str{sc}='c=zeros(ndec,1);';
    eval(str{sc});
    sc=sc+1;
    str{sc}='nc=decinfo(lmi,xg);';
    eval(str{sc});
    sc=sc+1;
    str{sc}='c(nc)=1;';
    eval(str{sc});
    sc=sc+1;
    options=E.options;
    str{sc}='[gain,xopt]=mincx(lmi,c,options);';
    sc=sc+1;
    str{sc}='gain=sqrt(gain);';

    if E.options(5)>100,
       save lmi_exe state* options c_f c_z

       if exist('struct1')
          save lmi_exe struct* -append
       end

       if exist('l_1_1')
          save lmi_exe l_* -append
       end

       if exist('r_1_1')
          save lmi_exe r_* -append
       end
   
       if exist('L_1_1')
          save lmi_exe L_* -append
       end

      if exist('R_1_1')
         save lmi_exe R_* -append
      end

      fid=fopen('lmi_exe.m','wt');
      for i=1:sc,
          fprintf(fid,[str{i} '\n']);
      end
      fclose(fid);
    end

    if vrb,
       disp(['  Solving with ' num2str(ndec) ' decision variables ...'])
    end
    
    eval(str{sc-1})
    eval(str{sc})

    global ABST
    ABST.xopt=xopt;
    ABST.E=E;
    ABST.c_f=c_f;
    ABST.c_z=c_z;

    if ~isempty(xopt),
       ABST.P=decs2mats(xopt,XP);
       ABST.xg=decs2mats(xopt,nc);
       % iqc_value;
    end

else
   export_date=date;
   str{1}=['%% Created by iqc_gain_tbx on ' export_date];
   str{2}='\n%% Load coefficient matrices ...';
   str{3}=['load lmi_exe'];

   % initialize YALMIP
   str{4}='\n%% Initialize YALMIP';
   str{5}='F = set([]);';   
   eval(str{5})

   % now define all variables
   if vrb,
      disp('  defining the original variables ...')
   end
   kvar=find(E.T==var);
   str{6}='\n%% Define multiplier variables ...';
   sc=7;              % string   counter
   vc=1;              % aux variable counter 
   vcs= num2str(vc);
   for k=1:length(kvar),
       ks   = num2str(k);
       varT = E.VT{kvar(k)}(1);
       varn = num2str(E.VT{kvar(k)}(2));
       varm = num2str(E.VT{kvar(k)}(3));
       switch varT
         case 1  % symmetric
           str{sc}=['X' ks '=sdpvar(' varn ',' varm ',''symmetric'');'];
           eval(str{sc});
           sc=sc+1;
         case 2  % rectangular
           str{sc}=['X' ks '=sdpvar(' varn ',' varm ',''full'');'];
           eval(str{sc});
           sc=sc+1;
         case 4  % skew symmetric
           str{sc}=['X' ks '=sdpvar(' varn ',' varm ',''skew'');'];
           eval(str{sc});
           sc=sc+1;
         case 3  % diagonal
           str{sc}=['X' ks '=zeros(' varn ',' varm ');'];
           eval(str{sc});
           sc = sc + 1; 
           for flcnt = 1: str2num(varn)
               str{sc}   = ['y' vcs '= sdpvar(1,1);'];
               str{sc+1} = ['X' ks '(' num2str(flcnt) ',' num2str(flcnt) ')=y' vcs ';'];
               eval(str{sc});
               eval(str{sc+1});
               sc = sc + 2; 
               vc = vc + 1;               
               vcs = num2str(vc);
           end
         case 5  % variable 
           str{sc}   = ['X' ks '=zeros(' varn ',' varm ');'];
           eval(str{sc});
           sc = sc + 1; 
           eval(['struct' ks '=E.L{' num2str(kvar(k)) '};']);

           M   = E.L{kvar(k)};
           Ma  = sort(M(:));
           Ms  = Ma(find(Ma>0));
           Msc = Ms(1);
           cnt = 1;
           for flcnt = 2: length(Ms)
               if Ms(flcnt) > Msc(cnt)           
                  Msc(cnt+1) = Ms(flcnt);
                  cnt = cnt + 1;
               end
           end

           for flcnt = 1: length(Msc)
               Mscvar    = num2str(Msc(flcnt));
               str{sc}   = ['y' vcs '= sdpvar(1,1);'];
               str{sc+1} = ['X' ks '=X' ks '+((abs(struct' ks ')==' Mscvar ').*sign(struct' ks '))*y' vcs ';']; 
               vc = vc + 1;               
               vcs = num2str(vc);
               eval(str{sc});
               eval(str{sc+1});
               sc = sc + 2; 
           end
       end
    end


    % now define all LMI's
    if vrb,
       disp('  defining the non-KYP LMIs ...')
    end
    kvar=find(E.T==lmi);
    str{sc}='\n%% Define non-KYP lmis ...';
    sc=sc+1;
    for k=1:length(kvar),
        kk=kvar(k);
        [nn,mm]=size(E.AB{kk});         % dimensions
        ks=num2str(k);
        nns=num2str(nn);
        vs=mm-nn;                       % LMI size
        vss=num2str(vs);
        str{sc}=['LMI' ks '=zeros(' num2str(mm) ');'];
        eval(str{sc});
        sc=sc+1;
        if nn>0,
           eval(['state' ks '=E.AB{kk};']);
           str{sc}=['P' ks '=sdpvar(' nns ',' nns ',''symmetric'');'];
           eval(str{sc});
           sc=sc+1;
           str{sc}=['LMI' ks '=LMI' ks '+[eye(' nns ');zeros(' vss ',' nns ')]*P' ...
                    ks '*state' ks '+state' ks '''*P' ks '*[eye(' nns ');zeros(' vss ',' ...
                    nns ')]'';'];
           eval(str{sc});
           sc=sc+1;
        end
        nh=0;
        ng=sum(E.X{kk}(3,:));
        for i=1:size(E.X{kk},2),
            is=num2str(i);
            nv=E.X{kk}(1,i);   % variable number
            if nv < 0
               nvs = num2str(-nv);
               istranspose = 1;
            else
               nvs = num2str(nv);
               istranspose = 0;
            end
            eval(['l_' ks '_' is '=E.C{kk}(ng+1:ng+E.X{kk}(2,i),:)'';']);
            eval(['r_' ks '_' is '=E.C{kk}(nh+1:nh+E.X{kk}(3,i),:);']);
            if istranspose
               str{sc}=['LMI' ks '=LMI' ks '+l_' ks '_' is '*X' nvs '''*' ...
                        'r_' ks '_' is '+r_' ks '_' is '''*X' nvs '*l_' ks '_' is ''';']; 
            else
               str{sc}=['LMI' ks '=LMI' ks '+l_' ks '_' is '*X' nvs '*' ...
                        'r_' ks '_' is '+r_' ks '_' is '''*X' nvs '''*l_' ks '_' is ''';'];  
            end
            eval(str{sc});
            sc=sc+1;
            nh=nh+E.X{kk}(3,i);
            ng=ng+E.X{kk}(2,i);
        end
        
 eval(['class(F)'])
 eval(['class(set(LMI' ks '<0)'])
 pause
        str{sc} = ['F = F + set(LMI' ks '<0);'];
        eval(str{sc});
        sc = sc + 1;
    end

    % now we define the "main" lmi 
    if vrb,
       disp('  defining the KYP LMIs ...')
    end
    str{sc}='\n%% KYP LMI terms...';
    sc=sc+1;
    nm=E.nlmi+1;         % LMI index  
    nms=num2str(nm);
    mm=E.ninputs;        % width of the 'B' matrix
    mms=num2str(mm);
    str{sc}=['KYPLMI=zeros(' num2str(E.nstates+E.ninputs) ');'];
    eval(str{sc});
    sc=sc+1;
    if E.nstates>0,
       nss=num2str(E.nstates);
       eval(['state0=E.ab;']);       
       str{sc}=['P=sdpvar(' nss ',' nss ',''symmetric'');'];
       eval(str{sc});
       sc=sc+1;
       str{sc}=['KYPLMI=KYPLMI+[eye(' nss ');zeros(' mms ',' nss ')]*P' ...
                 '*state0+state0''*P*[eye(' nss ');zeros(' mms ',' nss ')]'';'];
       eval(str{sc});
       sc=sc+1;   
    end
    for j=1:E.nsimple,
        js=num2str(j);
        nh=0;           % C-position counter for the right (C) terms 
        ng=0;           % C-position counter for the left (B) terms
        for i=1:size(E.X{E.simples(j)},2),
            is=num2str(i);
            vr = E.X{E.simples(j)}(1,i);
            if vr < 0
               vrs = num2str(-vr);
               istranspose = 1;
            else
               vrs = num2str(vr);
               istranspose = 0;
            end
            eval(['L_' js '_' is ...
                  '=E.B{E.simples(j)}(ng+1:ng+E.X{E.simples(j)}(2,i),:)'';'])
            eval(['R_' js '_' is ...
                  '=E.gqfm(j)*E.C{E.simples(j)}(nh+1:nh+' ...
                  'E.X{E.simples(j)}(3,i),:);'])
            if istranspose
               str{sc}=['KYPLMI=KYPLMI+L_' js '_' is '*X' vrs '''*R_' js '_' is ...
                        '+R_' js '_' is '''*X' vrs '*L_' js '_' is ''';'];          
            else
               str{sc}=['KYPLMI=KYPLMI+L_' js '_' is '*X' vrs '*R_' js '_' is ...
                        '+R_' js '_' is '''*X' vrs '''*L_' js '_' is ''';'];
            end
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
    str{sc}=['KYPLMI=KYPLMI+c_z''*c_z;'];
    eval(str{sc});
    sc=sc+1;
    str{sc}=['gamma_sq=sdpvar(1,1);'];
    eval(str{sc});
    sc=sc+1;
    str{sc}=['KYPLMI=KYPLMI-c_f''*gamma_sq*c_f;'];
    eval(str{sc});
    sc=sc+1;
    str{sc} = ['F = F + set(KYPLMI<0);'];
    eval(str{sc});
    sc = sc + 1;


    % solving the system of LMI's 
    str{sc}='\n%% Solving the system of LMIs ...';
    sc=sc+1;
    str{sc}=['ops = sdpsettings(''solver'',''' solver_option ''',''verbose'',0);'];
    eval(str{sc});
    sc=sc+1;
    str{sc}='sol  = solvesdp(F,gamma_sq,ops);';
    sc=sc+1;
    str{sc}='if sol.problem == 1 gain=[]; else gain = sqrt(double(gamma_sq)); end';
    sc=sc+1;
    str{sc}='varargout{1} = sol.solvertime;';

    if E.options(5)>100,
       save lmi_exe state* solver_option c_f c_z

       if exist('struct1')
          save lmi_exe struct* -append
       end

       if exist('l_1_1')
          save lmi_exe l_* -append
       end

       if exist('r_1_1')
          save lmi_exe r_* -append
       end
   
       if exist('L_1_1')
          save lmi_exe L_* -append
       end

      if exist('R_1_1')
         save lmi_exe R_* -append
      end

      fid=fopen('lmi_exe.m','wt');
      for i=1:sc,
          fprintf(fid,[str{i} '\n']);
      end
      fclose(fid);
    end

    %if vrb,
    %   disp(['  Solving with ' num2str(ndec) ' decision variables ...'])
    %end

    eval(str{sc-2})
    eval(str{sc-1})
    eval(str{sc})
    
    global ABST
    ABST.xopt=[];
    ABST.E=E;
    ABST.c_f=c_f;
    ABST.c_z=c_z;
end

% -------------------------------------------------
% the following function is for reducing the system
% program done by C. Kao on Jun. 15 1998
% last modified on June 16 1998
% last modified on June 17 1998
% last modified on June 21 1998
% --------------------------------------------------
 
function [E]=iqc_reduce(E)
% function [E]=iqc_reduce(E)
% 
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

    if A.log(a1,1)==inp & A.log(a2,1)==inp,
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
            c_index(forlp_counter);
            BigC=[BigC;E.C{c_index(forlp_counter)}];
            c_index_new=[c_index_new;c_index(forlp_counter)];
        end
    end
end
[mBigC,nBigC]=size(BigC);

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

    pos1=[E.POS(inp_replaced):E.POS(inp_replaced)+vs-1];
    pos2=[[ns+1:E.POS(inp_replaced)-1],[E.POS(inp_replaced)+vs:ns+ni]];
    elim=[elim,pos1];

    % check if this link has been processed before
    chk1=CD_lnk([cnt+1:cnt+vs],:);
    chk2=(chk1==zeros(size(chk1,1),size(chk1,2)));
    if ~all(chk2(:))
        D_replaced_lnk=CD_lnk([cnt+1:cnt+vs],pos1);
        r=rank(D_replaced_lnk,1e-6);
        if r~=vs,
           error([err_msg1,' ,non-invertable matrix found!'])
        else
           invD_replaced_lnk=inv(D_replaced_lnk);
           D_lnk=invD_replaced_lnk*(-1)*CD_lnk([cnt+1:cnt+vs],pos2);
           C_lnk=invD_replaced_lnk*(-1)*CD_lnk([cnt+1:cnt+vs],[1:ns]);
        end

% re-compute state space matrix 
% -----------------------------
        B_inp1=AB(:,pos1);
        B_reduced=AB(:,pos2)+B_inp1*D_lnk;
        A_reduced=AB(:,1:ns)+B_inp1*C_lnk;
        B_reduced_ext=zeros(mAB,nAB-ns);
        B_reduced_ext(:,[pos2-ns])=B_reduced;
        AB=[A_reduced,B_reduced_ext];
   
% re-compute output coefficient matrix
% ------------------------------------
        D_inp1=BigBC(:,pos1);
        D_reduced=BigBC(:,pos2)+D_inp1*D_lnk;
        C_reduced=BigBC(:,1:ns)+D_inp1*C_lnk;
        D_reduced_ext=zeros(mBigBC,nBigBC-ns);
        D_reduced_ext(:,[pos2-ns])=D_reduced;
        BigBC=[C_reduced,D_reduced_ext];

% re-compute the [C,D] matrix of 'Cx+Dw=0'
% ----------------------------------------
        D_inp1_CD_lnk=CD_lnk(:,pos1);
        D_reduced_CD_lnk=CD_lnk(:,pos2)+D_inp1_CD_lnk*D_lnk;
        C_reduced_CD_lnk=CD_lnk(:,[1:ns])+D_inp1_CD_lnk*C_lnk;
        D_reduced_ext_CD_lnk=zeros(mCD_lnk,nCD_lnk-ns);
        D_reduced_ext_CD_lnk(:,[pos2-ns])=D_reduced_CD_lnk;
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
   re_construct=[1:nBigC];
else
   re_construct=[];
   cnt=0;
   for forlp_counter=1:length(elim)
       if forlp_counter==length(elim)
          re_construct=[re_construct,[cnt+1:elim(forlp_counter)-1]];
          re_construct=[re_construct,[elim(forlp_counter)+1:nBigC]];
       else
           re_construct=[re_construct,[cnt+1:elim(forlp_counter)-1]];
           cnt=elim(forlp_counter);
       end
   end
end

% re-construct the state space matrix [A,B] and output matrices 
% in C field
% ---------------------------------------------------------------
E.ab=AB(:,re_construct);
BigB=BigBC([1:mBigB],re_construct);
BigC=BigBC([mBigB+1:mBigBC],re_construct);

cnt=0;
for forlp_counter=1:length(c_index_new)
    vs=size(E.C{c_index_new(forlp_counter)},1);
    E.C{c_index_new(forlp_counter)}=BigC([cnt+1:cnt+vs],:);
    cnt=cnt+vs;
end

cnt=0;
for forlp_counter=1:length(b_index_new)
    vs=size(E.B{b_index_new(forlp_counter)},1);
    E.B{b_index_new(forlp_counter)}=BigB([cnt+1:cnt+vs],:);
    cnt=cnt+vs;
end

cnt=0;
for forlp_counter=1:length(lnk_index)
    vs=size(E.C{lnk_index(forlp_counter)},1);
    E.C{lnk_index(forlp_counter)}=CD_lnk([cnt+1:cnt+vs],:);
    cnt=cnt+vs;
end

E.ninputs=E.ninputs-nreplaced;
