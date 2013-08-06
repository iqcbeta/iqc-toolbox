function y=num_op(x)
% function y=num_op(x)
%
% if x is a number, y is name of the operation #x
% if x is a string, y is the corresponding op. #
% if no such operation, error message
%
% Written by ameg@mit.edu,  last modified October 13, 1997

nop=30;          % number of supported operations

if ischar(x),
    switch x,
        case 'derivative', y=-3;
        case 'conversion', y=-1;
        case 'definition', y=0;
        case 'plus',       y=1;
        case 'minus',      y=2;
        case 'uminus',     y=3;
        case 'uplus',      y=4;
        case 'times',      y=5;
        case 'mtimes',     y=6;
        case 'rdivide',    y=7;
        case 'ldivide',    y=8;
        case 'mrdivide',   y=9;
        case 'mldivide',   y=10;
        case 'power',      y=11;
        case 'mpower',     y=12;
        case 'lt',         y=13;
        case 'gt',         y=14;
        case 'le',         y=15;
        case 'ge'     ,    y=16;
        case 'ne'    ,     y=17;
        case 'eq',         y=18;
        case 'and',        y=19;
        case 'or',         y=20;
        case 'not',        y=21;
        case 'colon',      y=22;
        case 'ctranspose', y=23;
        case 'transpose' , y=24;
        case 'display'   , y=25;
        case 'horzcat'   , y=26;
        case 'vertcat'   , y=27;
        case 'subsref'   , y=28;
        case 'subsasgn'  , y=29;
        case 'subsindex' , y=30;
        otherwise, error(['Unsupported operation *' x])
    end
else
    switch x,
        case -3,y='derivative';
        case -1,y='conversion';
        case 0, y='definition';
        case 1, y='plus';
        case 2, y='minus';
        case 3, y='uminus';
        case 4, y='uplus';
        case 5, y='times';
        case 6, y='mtimes';
        case 7, y='rdivide';
        case 8, y='ldivide';
        case 9, y='mrdivide';
        case 10,y='mldivide';
        case 11,y='power';
        case 12,y= 'mpower';
        case 13,y= 'lt';
        case 14,y= 'gt';
        case 15,y= 'le';
        case 16,y= 'ge';
        case 17,y= 'ne';
        case 18,y= 'eq';
        case 19,y= 'and';
        case 20,y= 'or';
        case 21,y= 'not';
        case 22,y= 'colon';
        case 23,y= 'ctranspose';
        case 24,y= 'transpose' ;
        case 25,y= 'display' ;
        case 26,y= 'horzcat';
        case 27,y= 'vertcat' ;
        case 28,y= 'subsref' ;
        case 29,y= 'subsasgn' ;
        case 30,y= 'subsindex' ;
        otherwise, error(['Unsupported operation #' num2str(x)])
    end
end