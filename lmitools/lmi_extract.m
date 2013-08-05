function E=lmi_extract(obj)
% function E=lmi_extract(obj)
%
% processes log file of "lmi" environment
% output is E - a structure
% which contains more detailed
% information about the system elements
%
% E.options :  same as ABST.options
% E.nlog    :  number of non-zero rows in ABST.log
% E.nlmivar :  total number of INDEPENDENT scalar decision variables
% E.nlmi    :  total number of LMI's
% E.nvar    :  total number of MATRIX decision variables
% E.nstates :  total number of states
% E.B{k}    :  lti system in a ss form (for "cst")
% E.L{k}    :  structure of a variable (for "var" only),
%              using the LMI Control Toolbox convention for "svar"
% E.G{k}    :  the "left" lti multiplier. For "cst", "var" and "lin"
%                   M=G'*X*H, where G,H are LTI, X is a block-diagonal
%                   variable.
% E.AB      :  local state-space models, used in frequency dependent LMI
% E.C{k}    :  output coefficients, used in LMI
% E.H{k}    :  the "right" lti multiplier for "var" and "lin"
% E.X{k}    :  the "block" structure of variables, a 3xNN array
%                  NN=# of blocks
%                  X(1,:) - matrix variable # (0 if constant, negative if transposed)
%                  X(2,:) - vertical size
%                  X(3,:) - horizontal size
%                  used in "cst", "var", "lin"
% E.T:      :  the first column of ABST.log (types of variables)
%
% to suppress the text output, use lmitbx_options to
%    set ABST.options(5)~=0 or 777
%
% written by cykao@mit.edu, last modified: Feb. 09 1999
%                           last modified: Apr. 29 1999
%                           last modified: Jul. 02 1999


% safety checks
global ABST
A=ABST;
if ~isfield(A,'name'),
   error('"abst" environment not defined')
end
if ~strcmp(A.name,'lmi'),
   error('This is not an "lmi" environment')
end

% vrb=0 to suppress the output
vrb=(A.options(5)==0)|(A.options(5)==777);

if nargin==0,
   obj=[];
   prob_feas=1;
elseif ~(A.log(double(obj),1)==2 | A.log(double(obj),1)==3)
   error('Bad cost function, not a ''var'' or a ''lin''')
elseif ~all(size(obj)==[1,1])
   error('cost function need to have size 1x1')
else prob_feas=0; end

% symbolic names for interior types:
cst=1;
var=2;
lin=3;
lmi=4;

if vrb,
   disp('lmi_extract: processing the abst log information ... ')
end

% first run: go through to process cst,var and lin entries
% and build the extra cell arrays 
%            (k is a log number)
%     B{k}  is an lti in the ss format (for cst)
%     L{k} - for variables
%     G{k} - left lti in ss format
%                 (nonempty for var,lin)
%     H{k} - right lti (for var,lin)
%     X{k} - list of variables (1st row - variable number from the tbx
%                               2nd and 3rd - vsize and hsize)

if vrb,
   disp('  Processing cst, var, and lin, counting ...')
end

AB{A.nlog}=[];
Q{1}=[];
L{1}=[];
G{1}=[];
H{1}=[];
X{1}=[];
C{1}=[];

nlmivar=0;        % number of LMI variables defined so far
nlmis=0;          % number of LMI's defined so far
nvar=0;
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
% --------------- type : constant ---------------------- %  
      case cst,                            % constant type
        switch op                          % what operation?
          case -1,                         % conversion
            B{k}=A.dat{a1};
          case 0,                          % definition
            switch df                      % definition of what ?
              case 7,                      % unit matrix
                B{k}=eye(vs);
              case 6,                      % zero matrix
                B{k}=zeros(vs,hs);
              otherwise
                error(['Invalid log entry ' num2str(k)])
            end
          case num_op('plus'),             % plus
            B{k}=B{a1}+B{a2};
          case num_op('minus'),            % minus
            B{k}=B{a1}-B{a2};
          case num_op('mtimes'),           % times
            B{k}=B{a1}*B{a2};  
          case num_op('uminus'),           % uminus
            B{k}=-B{a1};         
          case num_op('uplus'),            % uplus
            B{k}=B{a1};
          case num_op('vertcat'),          % vertcat
            B{k}=[B{a1};B{a2}];
          case num_op('horzcat'),          % horzcat
            B{k}=[B{a1},B{a2}];
          case num_op('ctranspose'),       % ctranspose
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
            B{k}=zeros(vs,hs);
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
        G{k}=eye(vs);
        H{k}=B{k};
        X{k}=[k*sqrt(-1);size(G{k},2);size(H{k},1)];
% --------------- type : variable ----------------------- % 
      case var,                             % variable type
        switch op                           % what operation?
          case 0,                           % definition
            switch df                       % definition of what ?
              case 1,                       % full symmetric
                svar=zeros(vs);             % structure matrix
                for i=1:vs,
                    for j=1:i,
                        nlmivar=nlmivar+1;
                        svar(j,i)=nlmivar;
                        svar(i,j)=nlmivar;
                    end
                end
              case 2,                       % full rectangular
                nlmivar=nlmivar+vs*hs;
                svar=reshape((nlmivar-vs*hs+1):nlmivar,hs,vs)';
              case 3,                       % diagonal
                nlmivar=nlmivar+vs;
                svar=diag((nlmivar-vs+1):nlmivar);
              case 4,                       % skew symmetric
                svar=zeros(vs);
                for i=1:vs,
                    for j=(i+1):vs,
                        nlmivar=nlmivar+1;
                        svar(i,j)=nlmivar;
                        svar(j,i)=-nlmivar;
                    end
                end
              case 5,                        % structured (defined by "variable")
                svar=ABST.dat{A.log(a1,6)};
                newvar=max(max(abs(svar)));
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
        G{k}=eye(vs);
        H{k}=eye(hs);
        X{k}=[nvar;vs;hs];
% --------------- type : linear combination ----------------------- %  
     case lin,                               % linear type
       switch op                             % produced by which operation ?
         case num_op('plus'),
           v1=size(G{a1},1);
           h1=size(H{a1},2);
           v2=size(G{a2},1);
           h2=size(H{a2},2);
           if (v1~=v2)|(h1~=h2),
              if (v1==1)&(h1==1),
                 G{k}=[ones(v2,1)*G{a1},G{a2}];  % horizontal concatenation
                 H{k}=[H{a1}*ones(1,h2);H{a2}];  % vertical concatenation
                 X{k}=[X{a1} X{a2}];
                 if v2==h2,
                    msg1='ambiguous scalar in matrix expression found in log # ';
                    msg2=[msg1,num2str(k),' : a+X here is intepreted as'];
                    msg3=[msg2,' a*ones(size(X))+X ! '];
                    warning(msg3);
                 end
              elseif (v2==1)&(h2==1),
                 G{k}=[G{a1},ones(v1,1)*G{a2}];  % horizontal concatenation
                 H{k}=[H{a1};H{a2}*ones(1,h1)];  % vertical concatenation
                 X{k}=[X{a1} X{a2}];
                 if v1==h1,
                    msg1='ambiguous scalar in matrix expression found in log # ';
                    msg2=[msg1,num2str(k),' : X+a here is intepreted as'];
                    msg3=[msg2,' X+a*ones(size(X)) ! '];
                    warning(msg3);
                 end
              else
                 error(['log # ' num2str(k) ': incompatible size in "plus"'])
              end
           else
              G{k}=[G{a1},G{a2}];  % horizontal concatenation
              H{k}=[H{a1};H{a2}];  % vertical concatenation
              X{k}=[X{a1} X{a2}];
           end
         case num_op('minus'),
           v1=size(G{a1},1);
           h1=size(H{a1},2);
           v2=size(G{a2},1);
           h2=size(H{a2},2);
           if (v1~=v2)|(h1~=h2),
              if (v1==1)&(h1==1),
                 G{k}=[ones(v2,1)*G{a1},-G{a2}];  % horizontal concatenation
                 H{k}=[H{a1}*ones(1,h2); H{a2}];  % vertical concatenation
                 X{k}=[X{a1} X{a2}];
                 if v2==h2,
                    msg1='ambiguous scalar in matrix expression found in log # ';
                    msg2=[msg1,num2str(k),' : a-X here is intepreted as'];
                    msg3=[msg2,' a*ones(size(X))-X ! '];
                    warning(msg3);
                 end
              elseif (v2==1)&(h2==1),
                 G{k}=[G{a1},-ones(v1,1)*G{a2}];  % horizontal concatenation
                 H{k}=[H{a1}; H{a2}*ones(1,h1)];  % vertical concatenation
                 X{k}=[X{a1} X{a2}];
                 if v1==h1,
                    msg1='ambiguous scalar in matrix expression found in log # ';
                    msg2=[msg1, num2str(k),' : X-a here is intepreted as'];
                    msg3=[msg2,' X-a*ones(size(X)) ! '];
                    warning(msg3);
                 end
              else
                 error(['log # ' num2str(k) ': incompatible size in "plus"'])
              end
           else
              G{k}=[G{a1},-G{a2}];  % horizontal concatenation
              H{k}=[H{a1}; H{a2}];  % vertical concatenation
              X{k}=[X{a1}  X{a2}];
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
                    X{k}=[X{a2}(1,1);h1;h1];
                    H{k}=eye(h1);
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
                    X{k}=[X{a1}(1,1);v2;v2];
                    G{k}=eye(v2);
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
           [v1,h1]=size(G{a1});
           [v2,h2]=size(G{a2});
           tmp1=[G{a1},zeros(v1,h2)];   % horizontal concatenation
           tmp2=[zeros(v2,h1),G{a2}];   % horizontal concatenation
           G{k}=[tmp1;tmp2];            % vertical concatenation
           X{k}=[X{a1} X{a2}];
           H{k}=[H{a1};H{a2}];          % vertical concatenation
         case num_op('horzcat'),
           [v1,h1]=size(H{a1});
           [v2,h2]=size(H{a2});
           tmp1=[H{a1},zeros(v1,h2)];   % horizontal concatenation
           tmp2=[zeros(v2,h1),H{a2}];   % horizontal concatenation
           H{k}=[tmp1;tmp2];            % vertical concatenation
           X{k}=[X{a1} X{a2}];
           G{k}=[G{a1},G{a2}];          % horizontal concatenation
         otherwise
           error(['Invalid lin log # ' num2str(k)])
       end
     end                             % finish of switch 'tp'
end                                  % finish first run loop

% second run: building the local state-space models
% storing state equations in a and output matrices in C{k}
% defining the lmi's lmis

nlmi=0;     % lmi counter
incnt=0;    % number of scalar inputs counted so far
stcnt=0;    % states counted so far
qfcnt=0;    % simple qfm counted so far
lmi_log=find(A.log(1:nlog,1)==lmi);

for k=1:length(lmi_log)
    tp=A.log(lmi_log(k),1);
    vs=A.log(lmi_log(k),2);
    hs=A.log(lmi_log(k),3);
    op=A.log(lmi_log(k),4);
    df=A.log(lmi_log(k),5);
    a1=A.log(lmi_log(k),6);
    a2=A.log(lmi_log(k),7);
    a3=A.log(lmi_log(k),8);
    if vs~=hs,
       error(['non-square lmi log #' num2str(lmi_log(k))])
    end
    s=[1 -1];
    if op==num_op('gt'),
       s=[-1 1];
    end
    % now building one big lmi GG*XX*HH<0

    v1=size(G{a1},1);
    h1=size(H{a1},2);
    v2=size(G{a2},1);
    h2=size(H{a2},2);
    if (v1~=v2)|(h1~=h2),
       if (v1==1)&(h1==1),
          done=0;
          if A.log(a1,1)==cst,
             if A.dat{A.log(a1,6)}==0,
                GG=s(2)*G{a2};
                HH=H{a2};
                X{lmi_log(k)}=X{a2};
                done=1;
             end
          end
          if ~done, 
             GG=[s(1)*ones(v2,1)*G{a1},s(2)*G{a2}];  % horizontal concatenation 
             HH=[H{a1}*ones(1,h2); H{a2}];           % vertical concatenation
             X{lmi_log(k)}=[X{a1}  X{a2}];
             if op==num_op('gt'),
                msg1='a > X';
                msg2='a*ones(size(X) > X';
             else
                msg1='a < X';
                msg2='a*ones(size(X) < X';
             end
             msg3='ambiguous scalar in matrix expression found in log # ';
             msg4=[msg3,num2str(lmi_log(k)),' : ', msg1,' here is intepreted as '];
             msg5=[msg4,msg2,'!'];
             warning(msg5);
          end

       elseif (v2==1)&(h2==1),
          done=0;
          if A.log(a2,1)==cst
             if A.dat{A.log(a2,6)}==0,
                GG=s(1)*G{a1};
                HH=H{a1};
                X{lmi_log(k)}=X{a1};
                done=1;
             end
          end
          if ~done,
             GG=[s(1)*G{a1},s(2)*ones(v1,1)*G{a2}];  % horizontal concatenation
             HH=[H{a1};H{a2}*ones(1,h1)];            % vertical concatenation
             X{lmi_log(k)}=[X{a1} X{a2}];
             if op==num_op('gt'),
                msg1='X > a';
                msg2='X > a*ones(size(X)';
             else
                msg1='X < a';
                msg2='X < a*ones(size(X)';
             end
             msg3='ambiguous scalar in matrix expression found in log # ';
             msg4=[msg3,num2str(lmi_log(k)),' : ',msg1,' here is intepreted as '];
             msg5=[msg4,msg2,'!'];
             warning(msg5);
          end
       end

    else
       GG=[s(1)*G{a1}, s(2)*G{a2}];  % horizontal concatenation
       HH=[H{a1}; H{a2}];            % vertical concatenation
       X{lmi_log(k)}=[X{a1}  X{a2}];
    end

    nh=size(HH,1); % # of outputs of HH
    DD=[HH;GG'];
    nlmis=nlmis+1;
    C{lmi_log(k)}=[DD];
    nlmi=nlmi+1;
    if vrb,
       disp(['    LMI #' sft(nlmi,5) 'size = ' sft(vs,5)])
    end
end

if prob_feas==0,      % whether optimization problem?
   k=double(obj);
   GG=G{k};
   HH=H{k};
   DD=[HH;GG'];
   C{k}=DD;
end
%%% end of second run %%%

% now define the output of fdlmi_extract
E.options=A.options;
E.nlog=nlog;
E.nlmivar=nlmivar;
E.nlmi=nlmi;
E.nvar=nvar;
E.nstate=stcnt;
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
   disp('  lmi_extract done OK')
   disp(' ')
end
