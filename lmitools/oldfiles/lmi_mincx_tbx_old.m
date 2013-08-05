function [copt,B,C,D,V,H,L]=lmi_mincx_tbx(obj)
% function lmi_mincx_tbx(obj)
%
% minimizes obj subject to the system of LMI's
% defined using the "abst" environment "lmi"
%
% obj must be class "abst" (var or lin)
% constant terms in obj are ignored in the case of lin
%
% if the problem is feasible,
% defines a global cell array ABSTSOLUTION of same size
% as ABST.log, containing the optimal doubles
%
%
% this function uses the LMI Control Toolbox by The MathWorks, Inc.
%
% Written by ameg@mit.edu,  last modified October 13, 1998
% Modified by cykao@mit.edu last modified Novenber 19, 1998


% first, we build a cell array, which,
% for each non-empty row of ABST.log, will contain the following:
%      constant entry -> the corresponding constant
%      variable or linear entry -> an Nx3 cell array 
%             { #of variable  left multiplier  right multiplier }
%      block or lmi - a rectangular matrix with reference numbers 
%      lmi - the same
%
% then we define a system of LMI's, remove unnecessary variables from
% memory, and run the LMI solver
%
% then for each entry in "log" we define the corresponding double,
% and produce the output

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

% check if the objective is a variable
%if A.log(double(varargin{1}),1)~=var,
%   error('First argument must be a variable')
%end
if ~isa(obj,'abst'),
   error('objective must be of class "abst"')
end

if ~ismember(A.log(double(obj)),[var lin]),
   error('objective must be of abstract class "var" or "lin"')
end

% check if the objective is a scalar
if ~isequal(size(obj),[1 1]),
   error('First argument must be a scalar')
end

setlmis([])       % initialize the LMI Control Toolbox
nlmivar=0;        % number of LMI variables defined so far
nlmis=0;          % number of LMI's defined so far

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
cfsd=0; % a counter for extra cell arrays in B

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
      C{k}={{eye(vs) nvar eye(hs)}};
      L{k}=svar;
   case lin,                     % linear type
      D{k}=k;
      V{k}=vs;
      H{k}=hs;
      switch op          % produced by which operation ?
      case num_op('plus'),
         B{k}=B{a1}+B{a2};
         if (A.log(a1,1)==var)&(A.log(a2,1)==var),
            B{k}=zeros(vs,hs);
         end
         C{k}={C{a1}{:} C{a2}{:}};
      case num_op('minus'),
         B{k}=B{a1}-B{a2};
         if (A.log(a1,1)==var)&(A.log(a2,1)==var),
            B{k}=zeros(vs,hs);
         end
         C{k}=C{a1};
         l0=length(C{a1});
         l1=length(C{a2});
         for l=1:l1,
            C{k}{l0+l}={C{a2}{l}{1} C{a2}{l}{2} -C{a2}{l}{3}};
         end
      case num_op('uplus'),
         B{k}=B{a1};
         if A.log(a1,1)==var,
            B{k}=zeros(vs,hs);
         end
         C{k}=C{a1};
      case num_op('uminus'),
         B{k}=-B{a1};
         if A.log(a1,1)==var,
            B{k}=zeros(vs,hs);
         end
         for l=1:length(C{a1}),
            C{k}{l}={C{a1}{l}{1} C{a1}{l}{2} -C{a1}{l}{3}};
         end
      case num_op('ctranspose'),
         B{k}=B{a1}';
         if A.log(a1,1)==var,
            B{k}=zeros(vs,hs);
         end
         for l=1:length(C{a1}),
            C{k}{l}={C{a1}{l}{3}' -C{a1}{l}{2} C{a1}{l}{1}'};
         end
      case num_op('mtimes'),
         B{k}=B{a1}*B{a2};
         if (A.log(a1,1)==var)|(A.log(a2,1)==var),
            B{k}=zeros(vs,hs);
         end
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
         elseif A.log(a1,1)==cst & H{a1}==sum(H{a2})
            column_pt=0;
            Dmatrix=[];
            for flcnt1=1:size(D{a2},2)
                B{A.nlog+cfsd+1}=B{a1}(:,column_pt+1:column_pt+H{a2}(flcnt1));
                Dmatrix=[Dmatrix,A.nlog+cfsd+1];
                column_pt=column_pt+H{a2}(flcnt1);
                cfsd=cfsd+1;
            end
            V{k}=[V{a1};V{a2}];
            D{k}=[Dmatrix;D{a2}];
            H{k}=H{a2};
         elseif A.log(a2,1)==cst & H{a2}==sum(H{a1})
            column_pt=0;
            Dmatrix=[];
            for flcnt1=1:size(D{a1},2)
                B{A.nlog+cfsd+1}=B{a2}(:,column_pt+1:column_pt+H{a1}(flcnt1));
                Dmatrix=[Dmatrix,A.nlog+cfsd+1];
                column_pt=column_pt+H{a1}(flcnt1);
                cfsd=cfsd+1;
            end            
            V{k}=[V{a1};V{a2}];
            D{k}=[D{a1};Dmatrix];
            H{k}=H{a1};
         else
            error(['lmi0: block sizes incompatible in ' num2str(k)])
         end
      case num_op('horzcat'),
         if isequal(V{a1},V{a2}),
            H{k}=[H{a1} H{a2}];
            D{k}=[D{a1} D{a2}];
            V{k}=V{a1};
         elseif A.log(a1,1)==cst & V{a1}==sum(V{a2})
            row_pt=0;
            Dmatrix=[];
            for flcnt1=1:size(D{a2},1)
                B{A.nlog+cfsd+1}=B{a1}(row_pt+1:row_pt+V{a2}(flcnt1),:);
                Dmatrix=[Dmatrix;A.nlog+cfsd+1];
                row_pt=row_pt+V{a2}(flcnt1);
                cfsd=cfsd+1;
            end
            H{k}=[H{a1},H{a2}];
            D{k}=[Dmatrix,D{a2}];
            V{k}=V{a2};
         elseif A.log(a2,1)==cst & V{a2}==sum(V{a1})
            row_pt=0;
            Dmatrix=[];
            for flcnt1=1:size(D{a1},1)
                B{A.nlog+cfsd+1}=B{a2}(row_pt+1:row_pt+V{a1}(flcnt1),:);
                Dmatrix=[Dmatrix;A.nlog+cfsd+1];
                row_pt=row_pt+V{a1}(flcnt1);
                cfsd=cfsd+1;
            end
            H{k}=[H{a1},H{a2}];
            D{k}=[D{a1},Dmatrix];
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
      if (~isequal(H{a1},H{a2},V{a1},V{a2}))&(~isequal(H{a1},V{a1},1))&(~isequal(H{a2},V{a2},1)),
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
        if size(D{p(s)},1)==size(D{p(s)},1) % check if block structure is square!
           for i=1:size(D{p(s)},1),       % vertical block number
               for j=1:size(D{p(s)},1),      % horizontal block number
                   d=D{p(s)}(i,j);
                   if i==j,
                       Bd=B{d}+B{d}';
                   else
                       Bd=B{d};
                   end
                   if ~isequal(Bd,zeros(A.log(d,2),A.log(d,3))),
                      lmiterm([q(s)*nlmis i j 0],Bd); % constant term
                   end
                   if d<=A.nlog % if d>A.nlog, there is only constant term!!
                      for l=1:length(C{d}),
                          if i~=j,
                             lmiterm([q(s)*nlmis i j C{d}{l}{2}],C{d}{l}{1},C{d}{l}{3});
                          else
                             lmiterm([q(s)*nlmis i j C{d}{l}{2}],C{d}{l}{1},C{d}{l}{3},'s');
                          end     % diagonal/non-diagonal terms
                      end      % loop over all terms in (i,j) block
                   end  
               end    % j-loop
           end        % i-loop
        end   
       end
      end
   otherwise
      error(['Invalid type in log ' num2str(k)])
   end                 
end

lmi=getlmis;
ndec=decnbr(lmi);
%c=zeros(ndec,1);
%nc=decinfo(lmi,C{double(varargin{1})}{1}{2});
%c(nc)=1;
c=zeros(ndec,1);
ko=double(obj);
for i=1:length(C{ko}),      % go through lmi terms
   vid=C{ko}{i}{2};         % "variable" part of the term
   cc=C{ko}{i}{1};          % term=cc*var*bb
   bb=C{ko}{i}{3};
   vn=abs(vid);             % variable number
   st=decinfo(lmi,vn);                % structure matrix
   if vid<0,                % if transposed ...
      st=st'; 
   end
   if isequal(size(st),[1 1]),
      c(abs(st))=c(abs(st))+cc*bb*sign(st);
   else
      for j=1:size(st,1),
         for l=1:size(st,2),
            c(abs(st(i,j)))=c(abs(st(i,j)))+cc(i)*bb(j)*sign(st(i,j));
         end
      end
   end
end


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
