function [w,Lambda]=iqc_popov_vect(v,sign)
% function [w,Lambda]=iqc_popov_vect(v,sign)
%
% Popov IQC for sector nonlinearities and parametric
% uncertainty (any size of v).
%
% w'*Lambda*(dv/dt)>0,
%
% where w=phi(v) in the case of nonlinearity and w=delta*v
% for the case with parametric uncertainty
%
% Here Lambda is symmetric nxn matrix with
%
% Lambda>0 if sign='+'
% Lambda<0 if sign='-'
% Lambda unconstrained if sign='0', which is default
%
% Work by Ulf Jonsson Nov 1997
% Last modified by cmj on 2013/5/4

if nargin<1
    disp_str(1)
end
if nargin<2
    sign='0';
end

global ABST

switch ABST.systemtype
    case 'continuous'
        n=size(v,1);
        u=derivative(v);
        w=signal(n);
        Lambda=symmetric(n);
        switch sign
            case '+'
                Lambda>0;
            case '-'
                Lambda<0;
            case '0'
            otherwise
                disp_str(69,'sign','''+'',''-'',''0'''),
        end
        w'*Lambda*u>0;
    case 'discrete'
        disp_str(70,'iqc_popov_vect','discrete')
end
