%%%%%%%%%%%%%%%%%%%%
%%                %%
%%  iqc_reduce    %%
%%                %%
%%%%%%%%%%%%%%%%%%%%

function [E]=iqc_reduce(E)
% function [E]=iqc_reduce(E)
% 
% internal function: 
% the purpose of this internal function is dealing with
% the 'link' stuff and reducing the system
% 
% written by cykao@mit.edu  June 22 1998

global ABST
A=ABST;

% symbolic names for interior types:
% ----------------------------------
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
% ------------------------------------------------------
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
% --------------------------------------------------------
b_index=qfm_index;

% define some useful constants
% ----------------------------
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
% -------------------------------------------------
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
% ------------------------------------------------------------------

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
% -------------------------------------------------------------
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
% -------------------------------------------------------------
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
