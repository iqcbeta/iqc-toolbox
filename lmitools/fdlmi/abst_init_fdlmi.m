function abst_init_fdlmi(type)
% function abst_init_fdlmi(type)
%
%   type: 0 or empty -> continuous system
%         1          -> discrete   system
%
% initializes the "abst" class environment for frequency dependent LMI handling
%
% internal classes:
%   1:  cst       (converted from "double","tf","ss","zpf", or "lti")
%   2:  var       (matrix  variable, can be declared with
%                  "symmetric","rectangular","diagonal"," skew", or "variable" )
%   3:  lin       (linear combinations of cst*var*cst)
%   4:  lmi       (inequalities with cst, var, lin)
%
% external classes:
%   double
%   tf
%   ss
%   zpk
%   lti
%
% Written by cykao@mit.edu, last modified: Feb 9, 1999
% Last modified by cmj on 2013/5/1

if nargin == 0
    type = 0;
end

% symbolic names for interior types:
cst=1;
var=2;
lin=3;
lmi=4;

% 保留舊 lmi 工具設定
global ABST ABSTSOLUTION
if isfield(ABST,'lmiparameter');
    lmitool=ABST.lmitool;
    lmiparameter=ABST.lmiparameter;
else
    chk=chklmitool;
    switch chk
        case {1,3}
            lmitool='lmilab';
            lmiparameter=[0 200 0 0 0];
        case 2
            lmitool='yalmip';
            lmiparameter=[];
    end
end

% 清除 ABST ABSTSOLUTION
ABST=[];
ABSTSOLUTION=[];

% 設定系統形式
switch type
    case 0
        ABST.systemtype = 'continuous';
    case 1
        ABST.systemtype = 'discrete';
    otherwise
        disp_str(69,'type','0,1')
end

ABST.name='fdlmi';              % environment name
ABST.mlog=8;                    % number of columns in ABST.log
ABST.log=zeros(100,ABST.mlog);  % will keep operations logic

% meaning of the columns of ABST.log:
% 1: interior class type (integer pointing to an element of ABST.cls)
% 2: vertical size
% 3: horizontal size
% 4: operation that led to forming the row (integer from "num_op")
% 5: comment - depends on the operation
%        definition:    definition flag (integer pointing to ABST.flag)
%        conversion:    which type converted from (points to ABST.ext)
% 6: first argument
%        conversion:    address of the data in ABST.dat
% 7: second argument
% 8: third argument: apparently, for subsassgn only
ABST.nlog=0;                    % number of USED rows in ABST.log
ABST.dat={};                    % will keep data
ABST.ndat=0;                    % no data yet
ABST.def={'symmetric'; ...      % for the use by "display"
    'rectangle'; ...
    'diagonal'; ...
    'skew'; ...
    'structure'; ...
    'zero'; ...
    'identity'; ...
    'input'; ...
    'vector'};
ABST.ext={'double'   1; ...
    'lti'      1; ...
    'ss'       1; ...
    'tf'       1; ...
    'zpk'      1};
ABST.next=5;                    % five external classes
ABST.cls={'constant'; ...
    'variable'; ...
    'linear'; ...
    'lmi'};
ABST.ncls=4;                    % number of interior types

% lmi toolbox
ABST.lmitool=lmitool;
ABST.lmiparameter=lmiparameter;

% tables for admissible operations follow
% rows correspond to the first argument type
% columns - to the second argument type
% entries correspond to the type of the result
% a zero entry means operation is not allowed

% unary operations available
ABST.uminus=zeros(1,ABST.ncls);
ABST.uplus=zeros(1,ABST.ncls);
ABST.ctranspose=zeros(1,ABST.ncls);
ABST.subsref=zeros(1,ABST.ncls);

% binary operations available
ABST.plus=zeros(ABST.ncls);
ABST.minus=zeros(ABST.ncls);
ABST.horzcat=zeros(ABST.ncls);
ABST.vertcat=zeros(ABST.ncls);
ABST.mtimes=zeros(ABST.ncls);
ABST.lt=zeros(ABST.ncls);
ABST.gt=zeros(ABST.ncls);
ABST.eq=zeros(ABST.ncls);
ABST.subsasgn=zeros(ABST.ncls);

% ************************************************************ cst
% admissible operations with "cst"
t=cst;
% unary operations
ABST.uplus(t)=t;
ABST.uminus(t)=t;
ABST.ctranspose(t)=t;
ABST.subsref(t)=t;
% binary operations
ABST.plus(t,[cst var lin])= ...
    [cst lin lin];
ABST.minus(t,[cst var lin])= ...
    [cst lin lin];
ABST.horzcat(t,[cst var lin])= ...
    [cst lin lin];
ABST.vertcat(t,[cst var lin])= ...
    [cst lin lin];
ABST.mtimes(t,[cst var lin])= ...
    [cst lin lin];
ABST.lt(t,[var lin])= ...
    [lmi lmi];
ABST.gt(t,[var lin])= ...
    [lmi lmi];
ABST.subsasgn(t,cst)= ...
    cst;

% ************************************************************** var
% admissible operations for "var"
t=var;
% unary operations
ABST.uplus(t)=lin;
ABST.uminus(t)=lin;
ABST.ctranspose(t)=lin;
ABST.subsref(t)=t;
% binary operations
ABST.plus(t,[cst var lin])= ...
    [lin lin lin];
ABST.minus(t,[cst var lin])= ...
    [lin lin lin];
ABST.vertcat(t,[cst var lin ])= ...
    [lin lin lin ];
ABST.horzcat(t,[cst var lin ])= ...
    [lin lin lin];
ABST.subsasgn(t,var)= ...
    var;
ABST.mtimes(t,cst)= ...
    lin;
ABST.lt(t,[cst,var,lin])= ...
    [lmi lmi lmi];
ABST.gt(t,[cst,var,lin])= ...
    [lmi lmi lmi];


% ************************************************************** lin
% admissible operations for "lin"
t=lin;
% unary operations
ABST.uplus(t)=lin;
ABST.uminus(t)=lin;
ABST.ctranspose(t)=lin;

% binary operations
ABST.plus(t,[cst var lin])= ...
    [lin lin lin];
ABST.minus(t,[cst var lin])= ...
    [lin lin lin];
ABST.vertcat(t,[cst var lin ])= ...
    [lin lin lin ];
ABST.horzcat(t,[cst var lin ])= ...
    [lin lin lin ];
ABST.mtimes(t,cst)= ...
    lin;
ABST.lt(t,[cst var lin])= ...
    [lmi lmi lmi];
ABST.gt(t,[cst var lin])= ...
    [lmi lmi lmi];

% ************************************************************** lmi
% no admissible operations for "lmi"

abst_const;
