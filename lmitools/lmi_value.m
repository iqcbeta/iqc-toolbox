function ABSTSOLUTION = lmi_value
% 
% this function, which can be run after an lmi-optimizer,
% such as "lmi_mincx_tbx", or "lmi_mincx_lp" calculates 
% "optimal" values of all "abst" elements and stores them 
% in the global variable ABSTSOLUTION (to be used later by 
% the "value" function)
%
% Written by cykao@mit.edu  last modified  Apr. 29 1999

global ABST
A=ABST;
clear global ABSTSOLUTION
global ABSTSOLUTION
ABSTSOLUTION={};

% safety check
if ~isfield(A,'name'),
   error('"abst" environment not defined')
end
if ~strcmp(A.name,'lmi'),
   error('This is not an "lmi" environment')
end
if ~isfield(A,'xopt'),
   error('You must run an optimizer (e.g. "lmi_mincx_tbx") first')
end

xopt=A.xopt;
E=A.E;

cst=1;
var=2;
lin=3;
lmi=4;
if ~isempty(xopt),
    X{1}=[];
    for flcnt=2:A.nlog,
        tp=A.log(flcnt,1);
        vs=A.log(flcnt,2);
        hs=A.log(flcnt,3);
        op=A.log(flcnt,4);
        df=A.log(flcnt,5);
        a1=A.log(flcnt,6);
        a2=A.log(flcnt,7);
        a3=A.log(flcnt,8);
        if tp==var,
           X{flcnt}=decs2mats(xopt,E.L{flcnt});
        else
           switch op
             case -1,                    % conversion
               X{flcnt}=A.dat{a1}; 
             case 0,                     % definition
               switch df                 % definition of what ?
                 case 7,                 % unit matrix
                   X{flcnt}=eye(vs);
                 case 6,                 % zero matrix
                   X{flcnt}=zeros(vs,hs);
                 otherwise
                   error(['Invalid constant definition log entry ' num2str(k)])
               end    
             case num_op('uplus'),
               X{flcnt}=X{a1};
             case num_op('uminus'),
               X{flcnt}=-X{a1};
             case num_op('ctranspose'),
               X{flcnt}=X{a1}';
             case num_op('plus'),        % plus
               X{flcnt}=X{a1}+X{a2};
             case num_op('minus'),       % minus
               X{flcnt}=X{a1}-X{a2};         
             case num_op('mtimes'),      % times
               X{flcnt}=X{a1}*X{a2};      
             case num_op('vertcat'),     % vertcat
               X{flcnt}=[X{a1};X{a2}]; 
             case num_op('horzcat'),     % horzcat
               X{flcnt}=[X{a1} X{a2}]; 
             case num_op('subsasgn'),
               X{flcnt}=X{a1};
               if df~=0,
                  X{flcnt}(a3,df)=X{a2};
               else
                  ind=A.dat{A.log(a3,6)};
                  for flcnt2=1:size(ind,1),
                      X{flcnt}(ind(flcnt2,1),ind(flcnt2,2))=X{a2}(ind(flcnt2,3),ind(flcnt2,4));
                  end
               end
             case num_op('subsref'),
               X{flcnt}=zeros(vs,hs);
               if a3~=0,
                  X{flcnt}=X{a1}(a2,a3);
               else
                  ind=A.dat{A.log(a2,6)};
                  for flcnt2=1:size(ind,1),
                      X{flcnt}(ind(flcnt2,1),ind(flcnt2,2))=X{a1}(ind(flcnt2,3),ind(flcnt2,4));
                  end
               end
             case num_op('lt'),
               X{flcnt}=(X{a1}<X{a2});
             case num_op('gt'),
               X{flcnt}=(X{a1}>X{a2});
             otherwise
               error(['Invalid log entry ' num2str(flcnt)])
             end
        end
    end
    ABSTSOLUTION=X;
end
