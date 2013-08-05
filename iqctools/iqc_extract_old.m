function E=iqc_extract(f,z)
% function E=iqc_extract(f,z)
%
% processes log file of iqc environment
% output is E - a structure 
% which contains more detailed
% information about the system elements
%
% E.options:  same as ABST.options
% E.nlog:     number of non-zero rows in ABST.log
% E.nlmivar:  total number of INDEPENDENT scalar decision variables
% E.nlmi:     total number of LMI's
% E.nvar:     total number of MATRIX decision variables
% E.ninputs:  total number of scalar inputs
% E.nstates:  total number of states (not counting the lmi-states)
% E.nsimple:  total number of simple quadratic forms
% E.ni:       input number in the log table
% E.no:       output number in the log table
% E.ab:       defines the state-space structure dx/dt=ab*[x;f]
%               where f is the vector of all inputs
% E.qfm:      coefficients of simple qf in global iqc (1-by-nsimple)
% E.simples:   indexes of simple qfm's
%
% E.B{k}:        lti system in a ss form (for "cst" and simple qfm)
% E.L{k}:        structure of a variable (for "var" only),
%                  using the LMI Control Toolbox convention for "svar"
% E.G{k}:        the "left" lti multiplier. For "var" and "lin"
%                  M=G'*X*H, where G,H are LTI, X is a
%                  block-diagonal variable. For "vsg" and "vcs"
%                  y(or z')=G'*X*CD*[x;f], where dx/dt=ab*[x;f] is the
%                  state-space model of the non-lmi dynamics
% E.AB:          local state-space models, used in lmi
% E.C{k}:        output coefficients, used in signals,lmi,simple qf
% E.H{k}:        the "right" lti multiplier for "var" and "lin"
% E.X{k}:        the "block" structure of variables, a 3xNN array
%                  NN=# of blocks
%                  x(1,:) - matrix variable # (negative if transposed)
%                  x(2,:) - vertical size
%                  x(3,:) - horizontal size
%                  used in var,lin, signals, and simple qforms
% E.Q{k}:        coefficients of simple qforms for "qfm" and "iqc"
%                  size 1-by-nsimple array
% E.T:           the first column of ABST.log (types of variables)
%                  
% to suppress the text output, use lmitbx_options to
%    set ABST.options(5)~=0 or 777
% Written by ameg@mit.edu, last modified November 10, 1997

% safety checks
global ABST
A=ABST;
if ~isfield(A,'name'), 
   error('"abst" environment not defined')
end
if ~strcmp(A.name,'iqc'),
   error('This is not an "iqc" environment')
end

% vrb=0 to suppress the output
vrb=(A.options(5)==0)|(A.options(5)==777);

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

if vrb,
   disp('iqc_extract: processing the abst log information ... ')
end

% first run: go through to count the number of inputs,
% states, and simple qfm's,
% resolve links between inputs
% process cst,var and lin entries
% building the extra cell arrays 
%            (k is a log number)
%     B{k}  is an lti in the ss format (for cst) 
%     L{k} - for variables
%     G{k} - left lti in ss format 
%                 (nonempty for var,lin,vsg,vcs)
%     H{k} - right lti (for var,lin,vsg,vcs)
%     X{k} - list of variables (1st row - variable number from the tbx
%                               2nd and 3rd - vsize and hsize)

if vrb,
   disp('  Processing cst, var, and lin, counting ...')
end

AB{1}=[];
Q{1}=[];
L{1}=[];
G{1}=[];
H{1}=[];
X{1}=[];
C{1}=[];

nlmivar=0;        % number of LMI variables defined so far
nlmis=0;          % number of LMI's defined so far
nvar=0;
ninputs=0;
nstates=0;
nsimple=0;
B{1}=[];          % k=1 corresponds to the empty abst element
G{1}=[];
H{1}=[];
X{1}=[];
L{1}=[];
for k=1:A.nlog,
   tp=A.log(k,1);
   vs=A.log(k,2);
   hs=A.log(k,3);
   op=A.log(k,4);
   df=A.log(k,5);
   a1=A.log(k,6);
   a2=A.log(k,7);
   a3=A.log(k,8);
   nlog=k;
   if tp==0,
      break
   end
   switch tp
   case inp,
      ninputs=ninputs+vs;     % count inputs
   case lnk,
      ninputs=ninputs-A.log(a1,2);    % reduce the count of inputs
      if (A.log(a1,2)~=A.log(a2,2)),
         error(['link incompatible at log #' num2str(k)])
      end
      if a1<a2,      % change younger log to indicate definition-by-link
         if A.log(a2,1)~=inp,
            error(['Younger in link not input at log #' num2str(k)])
         end
         A.log(a2,4)=-2;
         A.log(a2,6)=a1; 
      elseif a1>a2,
         if A.log(a1,1)~=inp,
            error(['Younger in link not input at log #' num2str(k)])
         end
         A.log(a1,4)=-2;
         A.log(a1,6)=a2;  
      else
         error(['link to itself detected at log # ' num2str(k)])
      end
   case qfm,
      if op==num_op('mtimes'),
         nsimple=nsimple+1;
      end
   case cst,                   % constant type
      switch op        % operation ?
      case -1,                    % conversion
         B{k}=ss(A.dat{a1});
      case 0,                     % definition
         switch df        % definition of what ?
         case 7,                      % unit matrix
            B{k}=ss(eye(vs));
         case 6,                      % zero matrix
            B{k}=ss(zeros(vs,hs));
         otherwise
            error(['Invalid log entry ' num2str(k)])
         end
      case num_op('plus'),        % plus
         B{k}=B{a1}+B{a2};
      case num_op('minus'),       % minus
         B{k}=B{a1}-B{a2};         
      case num_op('mtimes'),       % times
         B{k}=B{a1}*B{a2};         
      case num_op('uminus'),      % uminus
         B{k}=-B{a1};         
      case num_op('uplus'),        % uplus
         B{k}=B{a1};         
      case num_op('vertcat'),        % vertcat
         B{k}=mvrtct(B{a1},B{a2});         
      case num_op('horzcat'),        % horzcat
         B{k}=mhrzct(B{a1},B{a2}); 
      case num_op('ctranspose'),        % ctranspose
         B{k}=B{a1}';
      case num_op('subsasgn'),
         B{k}=B{a1};
         if df~=0,
            B{k}(a3,df)=B{a2};
         else
            ind=A.dat{A.log(a3,6)};
            for i=1:size(ind,1),
               B{k}(ind(i,1),ind(i,2))=B{a2}(ind(i,3),ind(i,4));
            end
         end
      case num_op('subsref'),
         B{k}=ss(zeros(vs,hs));
         if a3~=0,
            B{k}=B{a1}(a2,a3);
         else
            ind=A.dat{A.log(a2,6)};
            for i=1:size(ind,1),
               B{k}(ind(i,1),ind(i,2))=B{a1}(ind(i,3),ind(i,4));
            end
         end         
      otherwise
         error(['Invalid constant log entry ' num2str(k)])
      end
    case var,                   % variable type
      switch op
      case 0,
         switch df        % definition of what ?
         case 1,                      % full symmetric
            svar=zeros(vs);           % structure matrix
            for i=1:vs,
               for j=1:i,
                  nlmivar=nlmivar+1;
                  svar(j,i)=nlmivar;
                  svar(i,j)=nlmivar;
               end
            end
         case 2,                      % full rectangular
            nlmivar=nlmivar+vs*hs;
            svar=reshape((nlmivar-vs*hs+1):nlmivar,hs,vs)';
         case 3,                      % diagonal
            nlmivar=nlmivar+vs;
            svar=diag((nlmivar-vs+1):nlmivar);
         case 4,                      % skew symmetric
            svar=zeros(vs);
            for i=1:vs,
               for j=(i+1):vs,
                  nlmivar=nlmivar+1;
                  svar(i,j)=nlmivar;
                  svar(j,i)=-nlmivar;
               end
            end
         case 5,                      % structured (defined by "variable")
            svar=ABST.dat{A.log(a1,6)};
            newvar=max(abs(svar));
            svar=(nlmivar+abs(svar)).*sign(svar);
            nlmivar=nlmivar+newvar;
         otherwise
         error(['Invalid variable definition log entry ' num2str(k)])
         end
      case num_op('subsref'),
         svar=zeros(vs,hs);
         if a3~=0,
            svar=L{a1}(a2,a3);
         else
            ind=A.dat{A.log(a2,6)};
            for i=1:size(ind,1),
               svar(ind(i,1),ind(i,2))=L{a1}(ind(i,3),ind(i,4));
            end
         end
      case num_op('subsasgn'),
         svar=L{a1};
         if df~=0,
            svar(a3,df)=L{a2};
         else
            ind=A.dat{A.log(a3,6)};
            for i=1:size(ind,1),
               svar(ind(i,1),ind(i,2))=L{a2}(ind(i,3),ind(i,4));
            end
         end
      otherwise
         error(['Invalid variable log # ' num2str(k)])
      end  
      nvar=nvar+1;
      L{k}=svar; 
      G{k}=ss(eye(vs));
      H{k}=ss(eye(hs));
      X{k}=[nvar;vs;hs];
   case lin,                     % linear type
      switch op          % produced by which operation ?
      case num_op('plus'),
         if A.log(a1,1)==cst,
            G{k}=G{a2};
            H{k}=H{a2};
            X{k}=X{a2};
         elseif A.log(a2,1)==cst,
            G{k}=G{a1};
            H{k}=H{a1};
            X{k}=X{a1};
         else
            v1=size(G{a1},1);
            h1=size(H{a1},2);
            v2=size(G{a2},1);
            h2=size(H{a2},2);                 
            if (v1~=v2)|(h1~=h2),
               if (v1==1)&(h1==1),
                  G{k}=mhrzct(ones(v2,1)*G{a1},G{a2});
                  H{k}=mvrtct(H{a1}*ones(1,h2),H{a2});
               elseif (v2==1)&(h2==1),
                  G{k}=mhrzct(G{a1},ones(v1,1)*G{a2});
                  H{k}=mvrtct(H{a1},H{a2}*ones(1,h2));
               else   
                  error(['log # ' num2str(k) ': incompatible size in "plus"'])
               end   
            else
               G{k}=mhrzct(G{a1},G{a2});
               H{k}=mvrtct(H{a1},H{a2});
               X{k}=[X{a1} X{a2}];
            end
         end
      case num_op('minus'),
         if A.log(a1,1)==cst,
            G{k}=-G{a2};
            H{k}=H{a2};
            X{k}=X{a2};
         elseif A.log(a2,1)==cst,
            G{k}=G{a1};
            H{k}=H{a1};
            X{k}=X{a1};
         else
            v1=size(G{a1},1);
            h1=size(H{a1},2);
            v2=size(G{a2},1);
            h2=size(H{a2},2);                 
            if (v1~=v2)|(h1~=h2),
               if (v1==1)&(h1==1),
                  G{k}=mhrzct(ones(v2,1)*G{a1},-G{a2});
                  H{k}=mvrtct(H{a1}*ones(1,h2),H{a2});
               elseif (v2==1)&(h2==1),
                  G{k}=mhrzct(G{a1},-ones(v1,1)*G{a2});
                  H{k}=mhrzct(H{a1},H{a2}*ones(1,h2));
               else   
                  error(['log # ' num2str(k) ': incompatible size in "plus"'])
               end   
            else
               G{k}=mhrzct(G{a1},-G{a2});
               H{k}=mvrtct(H{a1},H{a2});
               X{k}=[X{a1} X{a2}];
            end  
         end
      case num_op('uplus'),
         G{k}=G{a1};
         H{k}=H{a1};
         X{k}=X{a1};
      case num_op('uminus'),
         G{k}=-G{a1};
         H{k}=H{a1};
         X{k}=X{a1}; 
      case num_op('ctranspose'), 
         G{k}=H{a1}';
         H{k}=G{a1}';
         X{k}=[-X{a1}(1,:);X{a1}(3,:);X{a1}(2,:)];
      case num_op('mtimes'),   
         % special care to multiplication by scalar variables  
         v1=A.log(a1,2);
         h1=A.log(a1,3);
         v2=A.log(a2,2);
         h2=A.log(a2,3);        
         if A.log(a1,1)==cst,
            if (h2==1)&(v2==1)&(h1~=v2),
               if (A.log(a2,1)==var),
                  X{k}=[X{a2}(2,1);h1;h1];
                  H{k}=tf(eye(h1));
                  G{k}=B{a1};
               else
                  error(['Invalid sizes in mtimes lin log # ' num2str(k)])
               end
            else
               X{k}=X{a2};
               G{k}=B{a1}*G{a2};
               H{k}=H{a2};
            end
         elseif A.log(a2,1)==cst,
            if (v2~=h1)&(h1==1)&(v1==1),
               if (A.log(a1,1)==var),
                  X{k}=[X{a1}(2,1);v2;v2];
                  G{k}=tf(eye(v2));
                  H{k}=B{a2};
               else
                  error(['Invalid sizes in mtimes lin log # ' num2str(k)])
               end
            else
               X{k}=X{a1};
               G{k}=G{a1};
               H{k}=H{a1}*B{a2};
            end          
         else
            error(['Invalid lin mtimes log entry ' num2str(k)])
         end
      case num_op('vertcat'),
         if A.log(a1,1)==cst,
            [v2,h2]=size(G{a2});
            G{k}=mvrtct(zeros(A.log(a1,2),h2),G{a2});
            X{k}=X{a2};
            H{k}=H{a2};
         elseif A.log(a2,1)==cst,
            [v1,h1]=size(G{a1});
            G{k}=mvrtct(G{a1},zeros(A.log(a2,2),h1));
            X{k}=X{a1};
            H{k}=H{a1};
         else   
            [v1,h1]=size(G{a1});
            [v2,h2]=size(G{a2});
            tmp1=mhrzct(G{a1},zeros(v1,h2));
            tmp2=mhrzct(zeros(v2,h1),G{a2});
            G{k}=mvrtct(tmp1,tmp2);
            X{k}=[X{a1} X{a2}];
            H{k}=mvrtct(H{a1},H{a2});
         end
      case num_op('horzcat'),
         if A.log(a1,1)==cst,
            [v2,h2]=size(H{a2});
            H{k}=mhrzct(zeros(v2,A.log(a1,3)),H{a2});
            X{k}=X{a2};
            G{k}=G{a2};
         elseif A.log(a2,1)==cst,
            [v1,h1]=size(H{a1});
            H{k}=mhrzct(H{a1},zeros(v1,A.log(a2,3)));
            X{k}=X{a1};
            G{k}=G{a1};
         else
            [v1,h1]=size(H{a1});
            [v2,h2]=size(H{a2});
            tmp1=mhrzct(H{a1},zeros(v1,h2));
            tmp2=mhrzct(zeros(v2,h1),H{a2});
            H{k}=mvrtct(tmp1,tmp2);
            X{k}=[X{a1} X{a2}];
            G{k}=mhrzct(G{a1},G{a2});
         end
      otherwise
         error(['Invalid lin log # ' num2str(k)])
      end      
   case sgn,
      if op==num_op('mtimes'),
         % disp(['k=' num2str(k) ';    sgn'])
         if (A.log(a1,2)==1)&(A.log(a1,3)==1),
            nstates=nstates+vs*size(ss(B{a1}),3);
         else
            nstates=nstates+size(ss(B{a1}),3);
         end
      end
   case vsg,
      if op==num_op('mtimes'),
         % disp(['k=' num2str(k) ';    vsg'])
         if A.log(a1,1)==cst,
            if (A.log(a1,2)==1)&(A.log(a1,3)==1),
               nstates=nstates+vs*size(ss(B{a1}),3);
            else
               nstates=nstates+size(ss(B{a1}),3);
            end
         else
            if (A.log(a1,2)==1)&(A.log(a1,3)==1),
               nstates=nstates+vs*(size(ss(G{a1}),3)+size(ss(H{a1}),3));
            else
               nstates=nstates+size(ss(G{a1}),3)+size(ss(H{a1}),3);
            end 
         end
      end
   case csg,
      if op==num_op('mtimes'),
         % disp(['k=' num2str(k) ';    csg'])
         if (A.log(a2,2)==1)&(A.log(a2,3)==1),
            nstates=nstates+hs*size(ss(B{a2}),3);
         else
            nstates=nstates+size(ss(B{a2}),3);
         end
      end      
   case vcs,
      if op==num_op('mtimes'),
       %  disp(['k=' num2str(k) ';    vcs'])
         if A.log(a2,1)==cst,
            if (A.log(a2,2)==1)&(A.log(a2,3)==1),         
               nstates=nstates+hs*size(ss(B{a2}),3);
            else
               nstates=nstates+size(ss(B{a2}),3);
            end
         else
            if (A.log(a2,2)==1)&(A.log(a2,3)==1),         
               nstates=nstates+hs*(size(ss(G{a2}),3)+size(ss(H{a2}),3));
            else
               nstates=nstates+size(ss(G{a2}),3)+size(ss(H{a2}),3);
            end 
         end
      end
   end   
end      % finish first run loop

if vrb,
   disp(['    scalar inputs:  ',num2str(ninputs)])
   disp(['    states:         ',num2str(nstates)])
   disp(['    simple q-forms: ',num2str(nsimple)])
   disp('  Processing signals and quadratic forms ...')
end
% second run: building the main state-space model
% storing state equations in a and output matrices in C{k}
% defining the lmi's lmis
% Q{k} - row vector of size nsimple - indicates signs of simple
%        q-forms in q-forms and iqc's

nlmi=0;     % lmi counter
incnt=0;    % number of scalar inputs counted so far
stcnt=0;    % states counted so far
qfcnt=0;    % simple qfm counted so far
a=zeros(nstates,nstates+ninputs);  
simples=zeros(1,nsimple);  % to keep addresses of simple qforms
gqfm=zeros(1,nsimple);     % coefficients of simple qfm in the qlobal iqc


for k=2:A.nlog,
   tp=A.log(k,1);
   vs=A.log(k,2);
   hs=A.log(k,3);
   op=A.log(k,4);
   df=A.log(k,5);
   a1=A.log(k,6);
   a2=A.log(k,7);
   a3=A.log(k,8);
   switch tp        % type ?
   case lmi,
      if vs~=hs,
         error(['non-square lmi log #' num2str(k)])
      end
      s=[1 -1];
      if op==num_op('gt'),
         s=[-1 1];
      end
      % now building one big lmi GG*XX*HH<0
      if A.log(a1,1)==cst,
         GG=s(2)*G{a2};
         HH=H{a2};
         X{k}=X{a2};
      elseif A.log(a2,1)==cst,
         GG=s(1)*G{a1};
         HH=H{a1};
         X{k}=X{a1};
      else
         if A.log(a1,2)~=A.log(a2,2),
            error('LMI with scalars are not implemented')
         else
            GG=mhrzct(s(1)*G{a1},s(2)*G{a2});
            HH=mvrtct(H{a1},H{a2});
         end
         X{k}=[X{a1} X{a2}];         
      end
      nh=size(HH,1); % # of outputs of HH 
      [AH,BH,CH,DH]=ssdata(HH);
      [AG,BG,CG,DG]=ssdata(GG);
      nG=size(AG,1);
      nH=size(AH,1);

      AA=[AH zeros(nH,nG);zeros(nG,nH) -AG'];
      BB=[BH;CG'];
      CC=[CH zeros(size(CH,1),nG);zeros(size(BG,2),nH) -BG'];
      DD=[DH;DG'];
      nlmis=nlmis+1;
      C{k}=[CC DD];
      AB{k}=[AA BB];
      nlmi=nlmi+1;
      if vrb,
      disp(['    LMI #' sft(nlmi,5) 'size = ' sft(vs,5) 'states: ' ...
            sft(size(CC,2),5) ]) 
      end
   case inp,
      switch op
      case 0,       % an original input
         C{k}=[zeros(vs,nstates+incnt) eye(vs) ...
               zeros(vs,ninputs-incnt-vs)] ;
         incnt=incnt+vs;
      case -2,      % a linked input
         C{k}=C{a1};
      otherwise
         error(['Invalid input log # ' num2str(k)])
      end
   case sgn,
      switch op
      case -3,                    % differentiation
         ddd=C{a1}(:,nstates+1:nstates+ninputs);
         ddd=ddd(:);
         if any(ddd~=0),
           error(['illegal differentiation log #' num2str(k)])
         end
         C{k}=C{a1}(:,1:nstates)*a;
      case num_op('plus'),        % plus
         if A.log(a1,2)==1,
            C{k}=ones(vs,1)*C{a1}+C{a2};
         elseif A.log(a2,2)==1,
            C{k}=C{a1}+ones(vs,1)*C{a2};
         else
            C{k}=C{a1}+C{a2};
         end
      case num_op('minus'),       % minus
         if A.log(a1,2)==1,
            C{k}=ones(vs,1)*C{a1}-C{a2};
         elseif A.log(a2,2)==1,
            C{k}=C{a1}-ones(vs,1)*C{a2};
         else
            C{k}=C{a1}-C{a2};
         end         
      case num_op('mtimes'),       % times
         if (A.log(a1,2)==1)&(A.log(a1,3)==1),
            ba1=B{a1}*eye(vs);
         elseif A.log(a2,2)==1,
            ba1=B{a1}*ones(size(B{a1},2),1);
         else
            ba1=B{a1};
         end
         [AA,BB,CC,DD]=ssdata(ba1);
         nn=size(AA,1);
         a(stcnt+1:stcnt+nn,stcnt+1:stcnt+nn)=AA;
         a(stcnt+1:stcnt+nn,:)=a(stcnt+1:stcnt+nn,:)+BB*C{a2};
         C{k}=DD*C{a2}; 
         C{k}(:,stcnt+1:stcnt+nn)=C{k}(:,stcnt+1:stcnt+nn)+CC;
         stcnt=stcnt+nn;
      case num_op('uminus'),      % uminus
         C{k}=-C{a1};         
      case num_op('uplus'),        % uplus
         C{k}=C{a1};         
      case num_op('vertcat'),        % vertcat
         if A.log(a1,1)==cst,
            C{k}=[zeros(A.log(a1,2),nstates+ninputs);C{a2}];
         elseif A.log(a2,1)==cst,
            C{k}=[C{a1};zeros(A.log(a2,2),nstates+ninputs)];
         else
            C{k}=[C{a1};C{a2}];
         end
      case num_op('ctranspose'),        % ctranspose
         C{k}=C{a1};        % no conjugation !!!
      case num_op('subsasgn'),
         C{k}=C{a1};
         if df~=0,               % single entry assignment
            C{k}(a3,:)=C{a2};
         else                    
            ind=A.dat{A.log(a3,6)};       % multiple entries assignment
            for i=1:size(ind,1),
               C{k}(ind(i,1),:)=C{a2}(ind(i,3),:);
            end
         end
      case num_op('subsref'),
         C{k}=zeros(vs,ninputs+nstates);
         if a3~=0,
            C{k}=C{a1}(a2,:);
         else
            ind=A.dat{A.log(a2,6)};
            for i=1:size(ind,1),
               C{k}(ind(i,1),:)=C{a1}(ind(i,3),:);
            end
         end         
      otherwise
         error(['Invalid signal log entry ' num2str(k)])
      end 
   case vsg,
      switch op          % produced by which operation ?
      case num_op('plus'),
         if A.log(a1,2)==1,
            G{k}=mhrzct(ones(vs,1)*G{a1},G{a2});
         elseif A.log(a2,2)==1,
            G{k}=mhrzct(G{a1},ones(vs,1)*G{a2});
         else
            G{k}=mhrzct(G{a1},G{a2});
         end
         C{k}=[C{a1};C{a2}];
         X{k}=[X{a1} X{a2}];
      case num_op('minus'),
         if A.log(a1,2)==1,
            G{k}=mhrzct(ones(vs,1)*G{a1},-G{a2});
         elseif A.log(a2,2)==1,
            G{k}=mhrzct(G{a1},-ones(vs,1)*G{a2});
         else
            G{k}=mhrzct(G{a1},-G{a2});
         end
         C{k}=[C{a1};C{a2}];
         X{k}=[X{a1} X{a2}];         
      case num_op('uplus'),
         C{k}=C{a1}; 
         G{k}=G{a1};
         X{k}=X{a1};
      case num_op('uminus'),
         C{k}=C{a1}; 
         G{k}=-G{a1};
         X{k}=X{a1};          
      case num_op('mtimes'),
         if A.log(a1,1)==cst,
            C{k}=C{a2};
            G{k}=B{a1}*G{a2};
            X{k}=X{a2};
         else
            sw=1; % switch to handle scalar variable exception (sw==0)
            if A.log(a2,2)==1,   % if the signal is scalar
               ca2=ones(A.log(a1,3),1)*C{a2};
            elseif A.log(a1,3)==A.log(a2,2), % if sizes are compatible
               ca2=C{a2};
            else
               if A.log(a1,1)==var, % scalar variable exception case
                  sw=0;
               else
                  error(['Incompatible size in vsg mtimes log # ' num2str(k)])
               end   
            end
            if sw,
               [AA,BB,CC,DD]=ssdata(H{a1});
               nn=size(AA,1);
               a(stcnt+1:stcnt+nn,stcnt+1:stcnt+nn)=AA;
               a(stcnt+1:stcnt+nn,:)=a(stcnt+1:stcnt+nn,:)+BB*ca2;
               C{k}=DD*ca2; 
               C{k}(:,stcnt+1:stcnt+nn)=C{k}(:,stcnt+1:stcnt+nn)+CC;
               stcnt=stcnt+nn;
               X{k}=X{a1};
               G{k}=G{a1};
            else
               nn=A.log(a2,2);
               C{k}=C{a2};
               G{k}=tf(eye(nn));
               X{k}=[X{a1}(1,:);nn;nn];
            end
         end
      otherwise
         error(['Invalid vsg log # ' num2str(k)])
      end
   case csg,
      switch op
      case num_op('plus'),        % plus
         if A.log(a1,3)==1,
            C{k}=ones(hs,1)*C{a1}+C{a2};
         elseif A.log(a2,3)==1,
            C{k}=C{a1}+ones(hs,1)*C{a2};
         else
            C{k}=C{a1}+C{a2};
         end         
      case num_op('minus'),       % minus
         if A.log(a1,3)==1,
            C{k}=ones(hs,1)*C{a1}-C{a2};
         elseif A.log(a2,3)==1,
            C{k}=C{a1}-ones(hs,1)*C{a2};
         else
            C{k}=C{a1}-C{a2};
         end         
      case num_op('mtimes'),       % times
         if (A.log(a2,2)==1)&(A.log(a2,3)==1),
            ba2=B{a2}'*eye(hs);
         elseif A.log(a1,3)==1,
            ba2=B{a2}'*ones(size(B{a2},1),1);
         else
            ba2=B{a2}';
         end
         [AA,BB,CC,DD]=ssdata(ba2);
         nn=size(AA,1);
         a(stcnt+1:stcnt+nn,stcnt+1:stcnt+nn)=AA;
         a(stcnt+1:stcnt+nn,:)=a(stcnt+1:stcnt+nn,:)+BB*C{a1};
         C{k}=DD*C{a1}; 
         C{k}(:,stcnt+1:stcnt+nn)=C{k}(:,stcnt+1:stcnt+nn)+CC;
         stcnt=stcnt+nn;
      case num_op('uminus'),      % uminus
         C{k}=-C{a1};         
      case num_op('uplus'),        % uplus
         C{k}=C{a1};         
      case num_op('horzcat'),        % still vertcat for C!!!
         if A.log(a1,1)==cst,
            C{k}=[zeros(A.log(a1,3),nstates+ninputs);C{a2}];
         elseif A.log(a2,1)==cst,
            C{k}=[C{a1};zeros(A.log(a2,3),nstates+ninputs)];
         else
            C{k}=[C{a1};C{a2}];
         end       
      case num_op('ctranspose'),        % ctranspose
         C{k}=C{a1};        % no conjugation !!!
      case num_op('subsasgn'),
         C{k}=C{a1};
         if df~=0,               % single entry assignment
            C{k}(a3,:)=C{a2};
         else                    
            ind=A.dat{A.log(a3,6)};       % multiple entries assignment
            for i=1:size(ind,1),
               C{k}(ind(i,1),:)=C{a2}(ind(i,3),:);
            end
         end
      case num_op('subsref'),
         C{k}=tf(zeros(vs,ninputs+nstates));
         if a3~=0,
            C{k}=C{a1}(a2,:);
         else
            ind=A.dat{A.log(a2,6)};
            for i=1:size(ind,1),
               C{k}(ind(i,1),:)=C{a1}(ind(i,3),:);
            end
         end         
      otherwise
         error(['Invalid co-signal log entry ' num2str(k)])
      end 
   case vcs,
      switch op          % produced by which operation ?
      case num_op('plus'),
         if A.log(a1,3)==1,
            G{k}=mhrzct(ones(hs,1)*G{a1},G{a2});
         elseif A.log(a2,3)==1,
            G{k}=mhrzct(G{a1},ones(hs,1)*G{a2});
         else
            G{k}=mhrzct(G{a1},G{a2});
         end
         C{k}=[C{a1};C{a2}];
         X{k}=[X{a1} X{a2}];
      case num_op('minus'),
         if A.log(a1,3)==1,
            G{k}=mhrzct(ones(hs,1)*G{a1},-G{a2});
         elseif A.log(a2,3)==1,
            G{k}=mhrzct(G{a1},-ones(hs,1)*G{a2});
         else
            G{k}=mhrzct(G{a1},-G{a2});
         end
         C{k}=[C{a1};C{a2}];
         X{k}=[X{a1} X{a2}];         
      case num_op('uplus'),
         C{k}=C{a1}; 
         G{k}=G{a1};
         X{k}=X{a1};
      case num_op('uminus'),
         C{k}=C{a1}; 
         G{k}=-G{a1};
         X{k}=X{a1}; 
      case num_op('ctranspose'),
         C{k}=C{a1};
         G{k}=G{a1};
         X{k}=X{a1};                 
      case num_op('mtimes'),
         if A.log(a2,1)==cst,
            C{k}=C{a1};
            G{k}=B{a2}'*G{a1};
            X{k}=X{a1};
         else
            sw=1; % switch to handle scalar variable exception
            if A.log(a1,3)==1,
               ca1=ones(A.log(a2,2),1)*C{a1};
            elseif A.log(a2,2)==A.log(a1,3),
               ca1=C{a1};
            else
               if A.log(a2,1)==var,
                  sw=0;
               else
                  error(['Incompatible size in vcs mtimes log # ' num2str(k)])
               end
            end
            if sw,
               [AA,BB,CC,DD]=ssdata(G{a2}');
               nn=size(AA,1);
               a(stcnt+1:stcnt+nn,stcnt+1:stcnt+nn)=AA;            
               a(stcnt+1:stcnt+nn,:)=a(stcnt+1:stcnt+nn,:)+BB*ca1;
               C{k}=DD*ca1; 
               C{k}(:,stcnt+1:stcnt+nn)=C{k}(:,stcnt+1:stcnt+nn)+CC;
               stcnt=stcnt+nn;
               X{k}=[-X{a2}(1,:);X{a2}(3,:);X{a2}(2,:)];
               G{k}=H{a2}';
            else
               nn=A.log(a1,3);
               C{k}=C{a1};
               G{k}=tf(eye(nn));
               X{k}=[-X{a2}(1,:);nn;nn];
            end
         end
      otherwise
         disp('oops')
         A.log(k,:)
         abst(k,0)
         error(['Invalid vcs log # ' num2str(k)])
      end
   case qfm,
      if (vs~=1)|(hs~=1),
         error(['non-scalar quadratic form log # ' num2str(k)])
      end
      switch op
      case num_op('mtimes'),
         if A.log(a1,1)==vcs,
            [AA,BB,CC,DD]=ssdata(G{a1}');
            nn=size(AA,1);
            a(stcnt+1:stcnt+nn,stcnt+1:stcnt+nn)=AA;            
            a(stcnt+1:stcnt+nn,:)=a(stcnt+1:stcnt+nn,:)+BB*C{a2};
            B{k}=DD*C{a2}; 
            B{k}(:,stcnt+1:stcnt+nn)=B{k}(:,stcnt+1:stcnt+nn)+CC;
            stcnt=stcnt+nn;
            X{k}=X{a1}; 
            C{k}=C{a1};
         elseif A.log(a2,1)==vsg,
            [AA,BB,CC,DD]=ssdata(G{a2}');
            nn=size(AA,1);
            a(stcnt+1:stcnt+nn,stcnt+1:stcnt+nn)=AA;            
            a(stcnt+1:stcnt+nn,:)=a(stcnt+1:stcnt+nn,:)+BB*C{a1};
            B{k}=DD*C{a1}; 
            B{k}(:,stcnt+1:stcnt+nn)=B{k}(:,stcnt+1:stcnt+nn)+CC;
            stcnt=stcnt+nn;
            X{k}=X{a2}; 
            C{k}=C{a2}; 
         else
            error(['Invalid qfm mtimes log # ' num2str(k)])
         end
         qfcnt=qfcnt+1;
         simples(qfcnt)=k;
         Q{k}=zeros(1,nsimple);
         Q{k}(qfcnt)=1;
      case num_op('uplus'),
         Q{k}=Q{a1};
      case num_op('uminus'),
         Q{k}=-Q{a1};
      case num_op('plus'),
         if A.log(a1,1)==cst,
            Q{k}=Q{a2};
         elseif A.log(a2,1)==cst,
            Q{k}=Q{a1};
         else
            Q{k}=Q{a1}+Q{a2};
         end
      case num_op('minus'),
         if A.log(a1,1)==cst,
            Q{k}=-Q{a2};
         elseif A.log(a2,1)==cst,
            Q{k}=Q{a1};
         else
            Q{k}=Q{a1}-Q{a2};
         end
      otherwise
         error(['Invalid q-form log # ' num2str(k)])
      end
   case iqc,  
      if A.log(a1,1)==cst,
         Q{k}=-Q{a2};
      elseif A.log(a2,1)==cst,
         Q{k}=Q{a1};
      else
         Q{k}=Q{a1}-Q{a2};
      end
      % positive sign corresponds to qfm>0
      if op==num_op('lt'),
         Q{k}=-Q{k};
      end
      gqfm=gqfm+Q{k};
   end
end          % end of second run

% check if the states count was correct
if stcnt~=nstates,
   error('the pre-calculated number of states is incorrect')
end


% now define the output of iqc_extract
E.options=A.options;
E.nlog=nlog;
E.nlmivar=nlmivar;
E.nlmi=nlmi;
E.nvar=nvar;
E.ninputs=ninputs;
E.nstates=nstates;
E.nsimple=nsimple;
E.ni=double(f);
E.no=double(z);
E.ab=a;
E.gqfm=gqfm;
E.simples=simples;
E.B=B;
E.L=L;
E.G=G;
E.AB=AB;
E.C=C;
E.H=H;
E.X=X;
E.Q=Q;
E.T=A.log(1:nlog,1);

if vrb,
   disp('  iqc_extract done OK')
   disp(' ')
end
