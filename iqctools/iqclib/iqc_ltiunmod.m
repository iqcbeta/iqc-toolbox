function [w,x]=iqc_ltiunmod(v,a,n,gain)
% Use: [w,x]=iqc_ltiunmod(v,a,n,gain)
%
% Defines w (size n) which relates to v by
%
%  w= Delta(s) v
%
% where Delta is an (uncertain) lti system
% with || Delta(s)||< gain.
%
% Arguments
%          a: (-) poles of multiplier
%          n:  output size
%          gain : norm bound
%
% Written by F. D'Amato,   last modified: June 11, 1998
% Last modified by cmj on 2013/5/4

ni=size(v,1);
if nargin < 4
    gain=1;
end
if nargin < 3
    n=ni;
end

global ABST

switch ABST.systemtype
    case 'continuous'
        if nargin < 2
            a=1;
        end
        k2=gain^2;
        w=signal(n);
        na=length(a);
        x=rectangular(na,1);
        x0=rectangular;
        sym_s=tf([1 0],1);
        p=x0;
        v'*k2*x0*v > w'*x0*w; %#ok<*VUNUS>
        for ndx=1:na;
            G=a(ndx)/(sym_s+a(ndx));
            vg=G*v;
            wg=G*w;
            p=p+x(ndx)*G;
            v'*k2*(x(ndx)*vg) > w'*(x(ndx)*wg);
        end
        p>0;
        x=[x0;x];
    case 'discrete'
        if nargin < 2
            a=-0.3679;
        end
        k2=gain^2;
        w=signal(n);
        na=length(a);
        x=rectangular(na,1);
        x0=rectangular;
        sym_z=tf([1 0],1,-1);
        p=x0;
        v'*k2*x0*v > w'*x0*w; %#ok<*VUNUS>
        for ndx=1:na;
            G=a(ndx)/(sym_z+a(ndx));
            vg=G*v;
            wg=G*w;
            p=p+x(ndx)*G;
            v'*k2*(x(ndx)*vg) > w'*(x(ndx)*wg);
        end
        p>0;
        x=[x0;x];
end
