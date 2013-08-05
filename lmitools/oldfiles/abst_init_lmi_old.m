function abst_init_lmi_old
% function abst_init_lmi_old 
% 
% initializes the "abst" class environment for lmi handling
% with a limited set of operations (subscripts and equalities
% not allowed, rectangular inequalities not allowed, block subdivisions
% must be "squarable", transposition of variables not allowed)
%
% though variables are real, constants can be complex-valued
%
% internal classes:
%   1:  cst       (converted from "double")
%   2:  var       (matrix variable, can be defined using
%                       symmetric(n)
%                       rectangular(n,m)
%                       diagonal(n)
%                       skew(n)
%                       variable(N)   etc.)
%   3:  lin       (linear combinations of variables and constants)
%   4:  lmi       (inequalities and equalities relating 2|3 and 1|2|3)
%   5:  blk       (block matrix function)
%
% external classes:
%   double
%
% Written by ameg@mit.edu,  last modified October 13, 1997

% symbolic names for interior types:
cst=1;
var=2;
lin=3;
lmi=4;
blk=5;

global ABST
if isfield(ABST,'options'),     % keep old options
   options=ABST.options;
else
   options=[0 200 0 0 0];
end

clear global ABST               % remove old trash
clear global ABSTSOLUTION
global ABST
global ABSTSOLUTION
ABSTSOLUTION={};
ABST.name='lmi';                % environment name
ABST.mlog=8;                    % number of columns in ABST.log
ABST.log=zeros(100,ABST.mlog);  % will keep operations logic

% meaning of the columns of ABST.log:
% 1: interior class type (integer pointing to an element of ABST.cls)
% 2: vertical size
% 3: horizontal size
% 4: operation that led to forming the row (integer from "num_op")
% 5: comment - depends on the operation
%        conversion:    which type converted from (points to ABST.ext)
% 6: first argument
%        conversion:    address of the data in ABST.dat
% 7: second argument
%

ABST.log(1,1:5)=[1 0 0 -1 1];  % reference for the empty object
ABST.nlog=1;                   % number of USED rows in ABST.log
ABST.dat={};                   % will keep data
ABST.ndat=0;                   % no data yet
ABST.ext={'double'   1}; 
ABST.next=1;                   % one external class 
ABST.def={'symmetric'; ...     % for the use by "display"
          'rectangle'; ...
          'diagonal'; ...
          'skew'; ...
          'structure'; ...
          'zero'; ...
          'identity'; ...
          'input'; ...
          'vector'};

ABST.cls={'cst'; ...
          'var'; ...
          'lin'; ...
          'lmi';
          'blk'};
% there will be 5 interior types
ABST.ncls=5;                 % 5 interior classes
ABST.options=options;

% tables for admissible operations follow
% rows correspond to the first argument type
% columns - to the second argument type
% entries correspond to the type of the result
% a zero entry means operation is not allowed

% unary operations available
%
%                     cst  var  lin  lmi  blk
% -------------------------------------------
ABST.uminus=       [  cst  lin  lin   0    0];
ABST.uplus=        [  cst  lin  lin   0    0];
ABST.ctranspose=   [  cst  lin  lin   0    0];
ABST.subsref=      [  cst  var   0    0    0];

% binary operations available
%
%                     cst  var  lin  lmi  blk
% -------------------------------------------
ABST.plus=         [  cst  lin  lin   0    0; ...  % cst
                      lin  lin  lin   0    0; ...  % var
                      lin  lin  lin   0    0; ...  % lin
                       0    0    0    0    0; ...  % lmi
                       0    0    0    0    0];     % blk

%                     cst  var  lin  lmi  blk
% -------------------------------------------
ABST.minus=        [  cst  lin  lin   0    0; ...  % cst
                      lin  lin  lin   0    0; ...  % var
                      lin  lin  lin   0    0; ...  % lin
                       0    0    0    0    0; ...  % lmi
                       0    0    0    0    0];     % blk

%                     cst  var  lin  lmi  blk
% -------------------------------------------
ABST.vertcat=      [  cst  blk  blk   0   blk; ...  % cst
                      blk  blk  blk   0   blk; ...  % var
                      blk  blk  blk   0   blk; ...  % lin
                       0    0    0    0    0 ; ...  % lmi
                      blk  blk  blk   0   blk];     % blk

%                     cst  var  lin  lmi  blk
% -------------------------------------------
ABST.horzcat=      [  cst  blk  blk   0   blk; ...  % cst
                      blk  blk  blk   0   blk; ...  % var
                      blk  blk  blk   0   blk; ...  % lin
                       0    0    0    0    0 ; ...  % lmi
                      blk  blk  blk   0   blk];     % blk

%                     cst  var  lin  lmi  blk
% -------------------------------------------
ABST.mtimes=       [  cst  lin  lin   0    0 ; ... % cst
                      lin   0    0    0    0 ; ... % var
                      lin   0    0    0    0 ; ... % lin
                       0    0    0    0    0 ; ... % lmi
                       0    0    0    0    0 ];    % blk

%                     cst  var  lin  lmi  blk
% -------------------------------------------
ABST.lt=           [   0   lmi  lmi   0   lmi; ... % cst
                      lmi  lmi  lmi   0    0 ; ... % cst
                      lmi  lmi  lmi   0    0 ; ... % cst
                       0    0    0    0    0 ; ... % lmi
                      lmi   0    0    0   lmi];    % blk

%                     cst  var  lin  lmi  blk
% -------------------------------------------
ABST.gt=           [   0   lmi  lmi   0   lmi; ... % cst
                      lmi  lmi  lmi   0    0 ; ... % cst
                      lmi  lmi  lmi   0    0 ; ... % cst
                       0    0    0    0    0 ; ... % lmi
                      lmi   0    0    0   lmi];    % blk

ABST.subsasgn=zeros(ABST.ncls);
ABST.subsasgn(cst,cst)=cst;
ABST.subsasgn(var,var)=var;

% now it's time to define useful constants
abst_const;
