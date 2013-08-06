function [cost,xopt]=fdlmi_mincx_tbx(obj)
% function [cost,xopt]=lmi_mincx_tbx(obj)
%
%  This program solves the LMI problem:
%
%       Minimize    c'x   subject to    L(x)  <  0
%       where L(x) is a matrix function. which is understood
%       as v'(L(x)+L(x)')v < 0 for all v in R^n.
%
%  The objective c'x (obj) has to be defined using the "abst"
%  environment "lmi". If the problem is feasible, a global cell
%  array ABSTSOLUTION, of same size as ABST.log, is created to
%  contain the optimal values of the decision variables.
%
% Written by cmj on 2013/4/29

global ABST

if ~isfield(ABST,'log'),
    disp_str(12)
end

if ~strcmp(ABST.name,'fdlmi')
    disp_str(15,'fdlmi')
end

lmitool = ABST.lmitool;

if nargin==1
    switch lmitool
        case 'lmilab'
            [cost,xopt] = fdlmi_mincx_lmilab(obj);
        case 'yalmip'
            [cost,xopt] = fdlmi_mincx_yalmip(obj);
        otherwise
            disp_str(55)
    end
elseif nargin==0
    switch lmitool
        case 'lmilab'
            [cost,xopt] = fdlmi_mincx_lmilab;
        case 'yalmip'
            [cost,xopt] = fdlmi_mincx_yalmip;
        otherwise
            disp_str(55)
    end
end
