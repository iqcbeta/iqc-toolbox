function iqc_value
% function iqc_value
% 
% this function, which can be run after an iqc-optimizer,
% such as "iqc_gain_tbx",
% calculates "optimal" values of all "abst" elements 
% and stores them in the global variable ABSTSOLUTION
% (to be used later by the "value" function)
%
% Written by ameg@mit.edu

% safety checks
global ABST
A=ABST;
if ~isfield(A,'name'),
   error('"abst" environment not defined')
end
if ~strcmp(A.name,'iqc'),
   error('This is not an "iqc" environment')
end
if ~isfield(A,'xopt'),
   error('You must run an optimizer (e.g. "iqc_gain_tbx") first')
end

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

xopt=A.xopt;
E=A.E;
statesystem=ss(E.ab(:,1:E.nstates), ...
   E.ab(:,E.nstates+1:E.nstates+E.ninputs), ...
   [eye(E.nstates);zeros(E.ninputs,E.nstates)], ...
   [zeros(E.nstates,E.ninputs);eye(E.ninputs)]);
Y{E.nlog}=[];
Y{1}=[];
kvar=find(E.T==var);
for k=1:length(kvar),        % values of variables
   Y{kvar(k)}=decs2mats(xopt,E.L{kvar(k)});
end
kvar=find(E.T==inp);
for k=1:length(kvar),        % ss: basic inputs to inp,sgn, csg
   Y{kvar(k)}=E.C{kvar(k)}*statesystem;
end

% main loop: defining all values (the Y matrix)
for k=2:A.nlog,
    tp=A.log(k,1);
    vs=A.log(k,2);
    hs=A.log(k,3);
    op=A.log(k,4);
    df=A.log(k,5);
    a1=A.log(k,6);
    a2=A.log(k,7);
    a3=A.log(k,8);
    switch  tp
      case var,
        Y{k}=decs2mats(xopt,E.L{k});
      case inp,
        Y{k}=tf(E.C{k}*statesystem);
      otherwise, 
        switch op
          case num_op('derivative'),
            [a,b,c,d]=ssdata(Y{a1});
            Y{k}=tf(ss(a,b,c*a,c*b));
          case -1,                        % conversion
            Y{k}=tf(A.dat{a1}); 
          case 0,                         % definition
            switch df                    % definition of what ?
              case 7,                      % unit matrix
                Y{k}=tf(eye(vs));
              case 6,                      % zero matrix
                Y{k}=tf(zeros(vs,hs));
              otherwise
                error(['Invalid constant definition log entry ' num2str(k)])
            end
          case num_op('uplus'),
            Y{k}=Y{a1};
          case num_op('uminus'),
            Y{k}=-Y{a1};
          case num_op('ctranspose'),
            Y{k}=Y{a1}';
          case num_op('plus'),         % plus
            Y{k}=Y{a1}+Y{a2};
          case num_op('minus'),        % minus
            Y{k}=Y{a1}-Y{a2};
          case num_op('mtimes'),       % times
            Y{k}=Y{a1}*Y{a2};
          case num_op('vertcat'),      % vertcat  
            if (A.log(a1,1)==cst)& ...
               ismember(A.log(a2,1),[inp sgn vsg]),
               Y{k}=tf(mvrtct(ss(zeros(A.log(a1,2), ...
               E.ninputs)),Y{a2}));
            elseif (A.log(a2,1)==cst)& ...
               ismember(A.log(a1,1),[inp sgn vsg]),
               Y{k}=tf(mvrtct(Y{a1},ss(zeros(A.log(a2,2), ...
               E.ninputs))));
            else
               Y{k}=tf(mvrtct(Y{a1},Y{a2}));
            end
          case num_op('horzcat'),      % horzcat
            if (A.log(a1,1)==cst)& ...
               ismember(A.log(a2,1),[csg vcs]),
               Y{k}=tf(mhrzct(ss(zeros(E.ninputs,A.log(a1,3))), ...
               Y{a2}));
            elseif (A.log(a2,1)==cst)& ...
               ismember(A.log(a1,1),[csg vcs]),
               Y{k}=tf(mhrzct(Y{a1},ss(zeros(E.ninputs, ...
               A.log(a2,3)))));
            else
               Y{k}=tf(mhrzct(Y{a1},Y{a2}));
            end
          case num_op('subsasgn'),
            if ismember(tp,[cst lin]),
               Y{k}=tf(Y{a1});
               if df~=0,
                  Y{k}(a3,df)=Y{a2};
               else
                  ind=A.dat{A.log(a3,6)};
                  for i=1:size(ind,1),
                      Y{k}(ind(i,1),ind(i,2))=Y{a2}(ind(i,3),ind(i,4));
                  end
               end
            elseif ismember(tp,[sgn vsg]),
               Y{k}=tf(Y{a1});
               if df~=0,
                  Y{k}(a3,:)=Y{a2};
               else
                  ind=A.dat{A.log(a3,6)};
                  for i=1:size(ind,1),
                      Y{k}(ind(i,1),:)=Y{a2}(ind(i,3),:);
                  end
               end
            else
               Y{k}=tf(Y{a1});
               if df~=0,
                  Y{k}(:,df)=Y{a2};
               else
                  ind=A.dat{A.log(a3,6)};
                  for i=1:size(ind,1),
                      Y{k}(:,ind(i,2))=Y{a2}(:,ind(i,4));
                  end
               end
            end
          case num_op('subsref'),
            if ismember(tp,[cst lin]),
               Y{k}=tf(zeros(vs,hs));
               if a3~=0,
                  Y{k}=Y{a1}(a2,a3);
               else
                  ind=A.dat{A.log(a2,6)};
                  for i=1:size(ind,1),
                      Y{k}(ind(i,1),ind(i,2))=Y{a1}(ind(i,3),ind(i,4));
                  end
               end
            elseif ismember(tp,[sgn vsg]),
               Y{k}=tf(zeros(vs,E.ninputs));
               if a3~=0,
                  Y{k}=Y{a1}(a2,:);
               else
                  ind=A.dat{A.log(a2,6)};
                  for i=1:size(ind,1),
                      Y{k}(ind(i,1),:)=Y{a1}(ind(i,3),:);
                  end
               end
            else
               Y{k}=tf(zeros(E.ninputs,hs));
               if a3~=0,
                  Y{k}=Y{a1}(a3,:);
               else
                  ind=A.dat{A.log(a2,6)};
                  for i=1:size(ind,1),
                      Y{k}(:,ind(i,2))=Y{a1}(:,ind(i,4));
                  end
               end
            end
          otherwise
          Y{k}=[];
        end
    end
end

clear global ABSTSOLUTION
global ABSTSOLUTION
ABSTSOLUTION=Y;
