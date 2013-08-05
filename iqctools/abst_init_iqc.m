function abst_init_iqc
% function abst_init_iqc
% 
% initializes the "abst" class environment for iqc handling
% 
%
% internal classes:
%   1:  cst       (converted from "double" and "lti")
%   2:  var       (matrix  variable, can be declared with
%       "symmetric","rectangular","diagonal"," skew", or "variable" )
%   3:  lin       (linear combinations of cst*var*cst)
%   4:  lmi       (inequalities with cst, var, lin)
%   5:  inp       (vector signal, can be declared with "input")
%   6:  sgn       (constant LTI transformations of inputs)
%   7:  vsg       (variable signal=lin-type transformation of a sgn)
%   8:  csg       (conjugated signal)
%   9:  vcs       (conjugated variable signals)
%  10:  qfm       (qfm=csg*vsg or vcs*sgn)
%  11:  iqc       (inequalities between qfm)
%  12:  lnk       (link between signals)
%
% external classes:
%   double
%   lti (also tf,ss,zpk)
%
% Written by ameg@mit.edu,  last modified October 13, 1997
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

global ABST
if isfield(ABST,'options'),
   options=ABST.options;
else
   options=[0 200 0 0 0];
end

clear global ABST               % remove old trash
clear global ABSTSOLUTION
global ABST
global ABSTSOLUTION

ABSTSOLUTION={};

ABST.name='iqc';         % environment name

ABST.mlog=8;            % number of columns in ABST.log
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

ABST.nlog=0;                % number of USED rows in ABST.log
ABST.dat={};                % will keep data
ABST.ndat=0;                % no data yet

ABST.def={'symmetric'; ...     % for the use by "display"
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
          'zpk'      1; ...
}; 
ABST.next=5;                % five external classes

ABST.cls={'constant'; ...
          'variable'; ...
          'linear'; ...
          'lmi'; ...
          'input'; ...
          'signal'; ...
          'varsignal'; ...
          'cosignal'; ...
          'varcosig';
          'q-form'; ...
          'iqc'; ...
          'link'};
% number of interior types
ABST.ncls=12;
% LMI control toolbox options
ABST.options=options;

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
ABST.horzcat(t,[cst var lin csg vcs])= ...
               [cst lin lin csg vcs];
ABST.vertcat(t,[cst var lin inp sgn vsg])= ...
               [cst lin lin sgn sgn vsg];
ABST.mtimes(t,[cst var lin inp sgn vsg])= ...
              [cst lin lin sgn sgn vsg];
ABST.lt(t,[var lin  qfm])= ...
          [lmi lmi  iqc];
ABST.gt(t,[var lin qfm])= ...
          [lmi lmi iqc];
ABST.subsasgn(t,[cst])= ...
                [cst];

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
ABST.subsasgn(t,[var])= ...
                [var];
ABST.mtimes(t,[cst inp sgn])= ...
              [lin vsg vsg ];
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
ABST.mtimes(t,[cst inp sgn])= ...
              [lin vsg vsg];
ABST.lt(t,[cst var lin])= ...
          [lmi lmi lmi];
ABST.gt(t,[cst var lin])= ...
          [lmi lmi lmi];


% ************************************************************** lmi
% no admissible operations for "lmi"


% ************************************************************** inp
% admissible operations for "inp"
t=inp;
% unary operations
ABST.uplus(t)=sgn;
ABST.uminus(t)=sgn;
ABST.ctranspose(t)=csg;
ABST.subsref(t)=sgn;

% binary operations
ABST.plus(t,[inp sgn])= ...
            [sgn sgn];
ABST.minus(t,[inp sgn])= ...
             [sgn sgn];
ABST.vertcat(t,[cst inp sgn])= ...
               [sgn sgn sgn];
ABST.eq(t,[inp sgn])= ...
          [lnk lnk];
ABST.subsasgn(t,[inp sgn])= ...
                [sgn sgn];


% ************************************************************** sgn
% admissible operations for "sgn"
t=sgn;
% unary operations
ABST.uplus(t)=sgn;
ABST.uminus(t)=sgn;
ABST.ctranspose(t)=csg;
ABST.subsref(t)=sgn;
% binary operations
ABST.plus(t,[inp sgn])= ...
            [sgn sgn];
ABST.minus(t,[inp sgn])= ...
             [sgn sgn];
ABST.vertcat(t,[cst inp sgn])= ...
               [sgn sgn sgn];
ABST.subsasgn(t,[inp sgn])= ...
                [sgn sgn];
ABST.eq(t,[inp])= ...
          [lnk];


% ************************************************************** csg
% admissible operations for "csg"
t=csg;
% unary operations
ABST.uplus(t)=csg;
ABST.uminus(t)=csg;
ABST.ctranspose(t)=sgn;
ABST.subsref(t)=csg;
% binary operations
ABST.plus(t,[csg])= ...
            [csg];
ABST.minus(t,[csg])= ...
            [csg];
ABST.mtimes(t,[cst var lin vsg])= ...
              [csg vcs vcs qfm];
ABST.horzcat(t,[cst csg])= ...
               [csg csg];
ABST.subsasgn(t,[csg])= ...
                [csg];


% ************************************************************** vsg
% admissible operations for "vsg"
t=vsg;
% unary operations
ABST.uplus(t)=vsg;
ABST.uminus(t)=vsg;
ABST.ctranspose(t)=vcs;

% binary operations
ABST.plus(t,[vsg])= ...
            [vsg];
ABST.minus(t,[vsg])= ...
            [vsg];

% ************************************************************** vcs
% admissible operations for "vcs"
t=vcs;
% unary operations
ABST.uplus(t)=vcs;
ABST.uminus(t)=vcs;
ABST.ctranspose(t)=vsg;

% binary operations
ABST.plus(t,[vcs])= ...
            [vcs];
ABST.minus(t,[vcs])= ...
             [vcs];
ABST.mtimes(t,[cst inp sgn])= ...
              [vcs qfm qfm];


% ************************************************************** iqc
% no admissible operations for "iqc"


% ************************************************************ qfm
% admissible operations with "qfm"
t=qfm;
% unary operations
ABST.uplus(t)=t;
ABST.uminus(t)=t;
% ABST.subsref(t)=t;
% binary operations
ABST.plus(t,[cst qfm])= ...
            [qfm qfm];
ABST.minus(t,[cst qfm])= ...
             [cst qfm];
ABST.gt(t,[cst qfm])= ...
          [iqc iqc];
ABST.lt(t,[cst qfm])= ...
          [iqc iqc];
ABST.eq(t,[cst qfm])= ...
          [iqc iqc];


% ************************************************************ link
% no admissible operations with "link"

abst_const;

