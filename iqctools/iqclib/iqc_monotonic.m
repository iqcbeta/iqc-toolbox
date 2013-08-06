function [w,xa,xb,xc,xd,xpop]=iqc_monotonic(v,a,op,k)
% function [w,xa,xb,xc,xd,xpop]=iqc_monotonic(v,a,op,k)
%
% defines iqc's for the relation
% w(t)=phi(v(t)),
% where phi is odd and its derivative stays in the interval [0,k].
%
% The IQC's have the form
%
% \int w'[1\pm 1/(1+s/a)](kv-w) > 0,                          ( Zames-Falb )
% \int (kv-w)'[1\pm 1/(1+s/a)]w > 0,
%
% and
% 0 < \dot{w}'*q*(v-w)                                        ( Popov )
%
% default N=1, a1=1, k=1
%
% op=1: shut off the Popov's IQC
% op=0; turn on  the Popov's IQC
% default: op=0
%
% Written by rantzer@control.lth.se, Sept 10, 1997. Last edited Nov 7, 1997.
% Modified by cykao@mit.edu Feb.2 1998, Jul. 12 1999
% Last modified by cmj on 2013/5/4

if nargin < 4
    k=1;
end
if nargin < 3
    op=0;
end

if nargin < 1
    disp_str(1)
end
if size(v,1)>1
    disp_str(8,'scalar')
end

global ABST

switch ABST.systemtype
    case 'continuous'
        if nargin < 2
            a=1;
        end
        s = tf([1 0],1);
        w = signal;
        r = k*v-w;
        m = length(a);
        xa = rectangular(m,1);
        xb = rectangular(m,1);
        xc = rectangular(m,1);
        xd = rectangular(m,1);
        for l=1:m,
            r1 = a(l)/(s+a(l))*r;
            w1 = a(l)/(s+a(l))*w;
            0 < w'*xa(l)*(r-r1); %#ok<*VUNUS>
            0 < w'*xb(l)*(r+r1);
            0 < r'*xc(l)*(w-w1);
            0 < r'*xd(l)*(w+w1);
            0 < xa(l);
            0 < xb(l);
            0 < xc(l);
            0 < xd(l);
        end,
        if ~op,
            u = derivative(v);
            xpop = symmetric;
            w'*xpop*u > 0;
        end
    case 'discrete'
        disp_str(70,'iqc_monotonic','discrete')
end
