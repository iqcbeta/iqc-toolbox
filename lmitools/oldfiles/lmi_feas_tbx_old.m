function lmi_feas_tbx
% function lmi_feas_tbx
%
% minimizes largest eigenvalue of the left part of LMI's
% defined using the "abst" environment "lmi"
%
%
% this function uses the LMI Control Toolbox by The MathWorks, Inc.
%
% Written by ameg@mit.edu,  last modified October 13, 1997

% safety checks
global ABST
A=ABST;
if ~isfield(A,'name'), 
   error('"abst" environment not defined')
end
if ~strcmp(A.name,'lmi'),
   error('This is not an "lmi" environment')
end

% symbolic names for interior types:
cst=1;
var=2;
lin=3;
lmi=4;
blk=5;


setlmis([])       % initialize the LMI Control Toolbox
nlmivar=0;        % number of LMI variables defined so far
nlmis=0;          % number of LMI's defined so far
[obj,nlmivar]=lmivar(1,[1 1]);


% building the extra cell arrays 
%     B  (for constants)
%     C  (for LMI terms)
%     D  (for block structure)
%     V  (vertical block sizes)
%     H  (horizontal block sizes)

% first element is an empty abstract element
B{1}=[];
C{1}={};
D{1}=0;
V{1}=0;
H{1}=0;
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
   case cst,                   % constant type
      C{k}={};
      D{k}=k;
      V{k}=vs;
      H{k}=hs;
      switch op        % operation ?
      case -1,                    % conversion
         B{k}=A.dat{a1}; 
      case 0,                     % definition
         switch df        % definition of what ?
         case 7,                      % unit matrix
            B{k}=eye(vs);
         case 6,                      % zero matrix
            B{k}=zeros(vs,hs);
         otherwise
            error(['Invalid constant definition log entry ' num2str(k)])
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
         B{k}=[B{a1};B{a2}];         
      case num_op('horzcat'),        % horzcat
         B{k}=[B{a1} B{a2}]; 
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
         B{k}=tf(zeros(vs,hs));
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
      B{k}=0;
      D{k}=k;
      V{k}=vs;
      H{k}=hs;
      switch op
         case 0,
         switch df        % definition of what ?
         case 1,                      % full symmetric
            [nvar,nlmivar,svar]=lmivar(1,[vs 1]);
         case 2,                      % full rectangular
            [nvar,nlmivar,svar]=lmivar(2,[vs hs]);
         case 3,                      % diagonal
            [nvar,nlmivar,svar]=lmivar(1,ones(vs,2));
         case 4,                      % skew symmetric
            struct=zeros(vs);
            for i=1:vs,
               for j=(i+1):vs,
                  nlmivar=nlmivar+1;
                  struct(i,j)=nlmivar;
                  struct(j,i)=-nlmivar;
               end
            end
            [nvar,nlmivar,svar]=lmivar(3,struct);
         case 5,
            struct=ABST.dat{a1};
            newvar=max(struct(:));
            struct=struct+nlmivar*(struct>0);
            [nvar,nlmivar,svar]=lmivar(3,struct);
         otherwise
         error(['Invalid variable definition log entry ' num2str(k)])
         end
      case num_op('subsref'),
         struct=zeros(vs,hs);
         if a3~=0,
            struct=L{a1}(a2,a3);
         else
            ind=A.dat{A.log(a2,6)};
            for i=1:size(ind,1),
               struct(ind(i,1),ind(i,2))=L{a1}(ind(i,3),ind(i,4));
            end
         end
         [nvar,nlmivar,svar]=lmivar(3,struct);
      case num_op('subsasgn'),
         struct=L{a1};
         if df~=0,
            struct(a3,df)=L{a2};
         else
            ind=A.dat{A.log(a3,6)};
            for i=1:size(ind,1),
               struct(ind(i,1),ind(i,2))=L{a2}(ind(i,3),ind(i,4));
            end
         end
         [nvar,nlmivar,svar]=lmivar(3,struct);
      otherwise
         error(['Invalid variable log # ' num2str(k)])
      end  
      C{k}={{1 nvar 1}};
      L{k}=svar;        
   case lin,                     % linear type
      D{k}=k;
      V{k}=vs;
      H{k}=hs;
      switch op          % produced by which operation ?
      case num_op('plus'),
         B{k}=B{a1}+B{a2};
         C{k}={C{a1}{:} C{a2}{:}};
      case num_op('minus'),
         B{k}=B{a1}-B{a2};
         C{k}=C{a1};
         l0=length(C{a1});
         l1=length(C{a2});
         for l=1:l1,
            C{k}{l0+l}={C{a2}{l}{1} C{a2}{l}{2} -C{a2}{l}{3}};
         end
      case num_op('uplus'),
         B{k}=B{a1};
         C{k}=C{a1}; 
      case num_op('uminus'),
         B{k}=-B{a1};
         for l=1:length(C{a1}),
            C{k}{l}={C{a1}{l}{1} C{a1}{l}{2} -C{a1}{l}{3}}; 
         end
      case num_op('ctranspose'), 
         B{k}=B{a1}';
         for l=1:length(C{a1}),
            C{k}{l}={C{a1}{l}{3}' -C{a1}{l}{2} C{a1}{l}{1}'}; 
         end             
      case num_op('mtimes'),
         B{k}=B{a1}*B{a2};
         if length(C{a2})==0,
            for l=1:length(C{a1}),
               C{k}{l}={C{a1}{l}{1} C{a1}{l}{2} C{a1}{l}{3}*B{a2}}; 
            end 
         elseif length(C{a1})==0,
            for l=1:length(C{a2}),
               C{k}{l}={B{a1}*C{a2}{l}{1} C{a2}{l}{2} C{a2}{l}{3}}; 
            end            
         else
            error(['Invalid multiplication entries ' num2str(k)])
         end
      otherwise
         error(['Invalid linear entry ' num2str(k)])
      end
   case blk,           % block matrix type
      switch op
      case num_op('vertcat'),
         if isequal(H{a1},H{a2}),
            V{k}=[V{a1};V{a2}];
            D{k}=[D{a1};D{a2}];
            H{k}=H{a1};
         else
            error(['lmi0: block sizes incompatible in ' num2str(k)])
         end
      case num_op('horzcat'),
         if isequal(V{a1},V{a2}),
            H{k}=[H{a1} H{a2}];
            D{k}=[D{a1} D{a2}];
            V{k}=V{a1};
         else
            error(['lmi0: block sizes incompatible in ' num2str(k)])
         end  
      otherwise
         error(['Invalid block entry ' num2str(k)]) 
      end
   case lmi,           % lmi type
      nlmis=nlmis+1;
      if num_op('lt')==op,
         p=[a1 a2];
      else
         p=[a2 a1];
      end
      q=[1 -1];
      if (~isequal(H{a1},H{a2},V{a1},V{a2}))&(~isequal(H{a1},V{a1},1)) ...
            &(~isequal(H{a2},V{a2},1)),
         error(['Incompatible LMI sides in ' num2str(k)])
      end
      
      for s=1:2,
       if (A.log(p(s),1)==cst)&(A.log(p(s),2)==1)&(A.log(p(s),3)==1), 
        for i=1:size(D{p(3-s)},1),       % vertical block number
         if B{p(s)}~=0,
          lmiterm([q(s)*nlmis i i 0],2*B{p(s)});
         end
        end
       else 
        % going through the terms
        for i=1:size(D{p(s)},1),       % vertical block number
         for j=1:size(D{p(s)},1),      % horizontal block number
          d=D{p(s)}(i,j);
          if i==j,
             Bd=B{d}+B{d}';
             if s==1,
                lmiterm([-nlmis i i obj],2,1);
             end
          else
            Bd=B{d};
          end
          if ~isequal(Bd,zeros(A.log(d,2),A.log(d,3))),
           lmiterm([q(s)*nlmis i j 0],Bd); % constant term
          end
          for l=1:length(C{d}),
           if i~=j,
            lmiterm([q(s)*nlmis i j C{d}{l}{2}],C{d}{l}{1},C{d}{l}{3});
           else
            lmiterm([q(s)*nlmis i j C{d}{l}{2}],C{d}{l}{1},C{d}{l}{3},'s');
           end     % diagonal/non-diagonal terms
          end      % loop over all terms in (i,j) block
         end       % j-loop
        end        % i-loop
       end
      end
   otherwise
      error(['Invalid type in log ' num2str(k)])
   end                 
end

lmi=getlmis;
ndec=decnbr(lmi);
c=zeros(ndec,1);
nc=decinfo(lmi,obj);
c(nc)=1;
[copt,xopt]=mincx(lmi,c,A.options);

clear global ABSTSOLUTION
global ABSTSOLUTION

ABSTSOLUTION={};
global ABSTSOLUTION
if ~isempty(xopt),
   X{1}=[];
   for k=2:A.nlog,
      tp=A.log(k,1);
      vs=A.log(k,2);
      hs=A.log(k,3);
      op=A.log(k,4);
      df=A.log(k,5);
      a1=A.log(k,6);
      a2=A.log(k,7);
      a3=A.log(k,8);
      if tp==var,
         X{k}=dec2mat(lmi,xopt,C{k}{1}{2});
      else
         switch op
         case -1,                    % conversion
            X{k}=A.dat{a1}; 
         case 0,                     % definition
            switch df        % definition of what ?
            case 7,                      % unit matrix
               X{k}=eye(vs);
            case 6,                      % zero matrix
               X{k}=zeros(vs,hs);
            otherwise
               error(['Invalid constant definition log entry ' num2str(k)])
            end    
         case num_op('uplus'),
            X{k}=X{a1};
         case num_op('uminus'),
            X{k}=-X{a1};
         case num_op('ctranspose'),
            X{k}=X{a1}';
         case num_op('plus'),        % plus
            X{k}=X{a1}+X{a2};
         case num_op('minus'),       % minus
            X{k}=X{a1}-X{a2};         
         case num_op('mtimes'),       % times
            X{k}=X{a1}*X{a2};               
         case num_op('vertcat'),        % vertcat
            X{k}=[X{a1};X{a2}];         
         case num_op('horzcat'),        % horzcat
            X{k}=[X{a1} X{a2}]; 
         case num_op('subsasgn'),
            X{k}=X{a1};
            if df~=0,
               X{k}(a3,df)=X{a2};
            else
               ind=A.dat{A.log(a3,6)};
               for i=1:size(ind,1),
                  X{k}(ind(i,1),ind(i,2))=X{a2}(ind(i,3),ind(i,4));
               end
            end
         case num_op('subsref'),
            X{k}=tf(zeros(vs,hs));
            if a3~=0,
               X{k}=X{a1}(a2,a3);
            else
               ind=A.dat{A.log(a2,6)};
               for i=1:size(ind,1),
                  X{k}(ind(i,1),ind(i,2))=X{a1}(ind(i,3),ind(i,4));
               end
            end  
         case num_op('lt'),
            X{k}=(X{a1}<X{a2});
         case num_op('gt'),
            X{k}=(X{a1}>X{a2});
         otherwise
            error(['Invalid log entry ' num2str(k)])
         end
      end
   end
   ABSTSOLUTION=X;
end

      

      
            
            
      
   